/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 9, 2015
 *      			Author: Minkwan Kim
 *
 * processing_grid.cpp
 * 			-  
 *  
 */

#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <limits>

//#include "../../Common/include/D_text.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../../MATRIX/include/math_misc.hpp"

#include "../include/constants_grid.hpp"
#include "../include/node_class.hpp"
#include "../include/face_class.hpp"
#include "../include/cell_class.hpp"
#include "../include/grid_class.hpp"



/*
 * ================================================================================
 * 		Adjust node data
 * ================================================================================
 */

void processing_node_data(unsigned int DIM, unsigned int	NNM, vector<NODE_CLASS *>	&nodes, double mesh_factor, bool is_axisymmetric)
{
	for (int n = 0; n <= NNM-1; n++)
	{
		for (int k = 0; k <= DIM-1; k++)
		{
			nodes[n]->x[k]	/= mesh_factor;
		}
	}


	if (is_axisymmetric	== true)
	{
		for (int n = 0; n <= NNM-1; n++)
		{
			nodes[n]->x[1]	+= 1.0e-9;
		}
	}
}



/*
 * ================================================================================
 * 		Process face data
 * ================================================================================
 */
void processing_face_data(unsigned int DIM, unsigned int NFM, vector<NODE_CLASS *>	&nodes, vector<int> &whereis_node, vector<FACE_CLASS *>	&faces)
{
	for (int f = 0; f <= NFM-1; f++)
	{
		// 1. Calculate Xc
		for (int k = 0; k <= DIM-1; k++)
		{
			faces[f]->x[k]	= 0.0;

			for (int n	= 0; n <= faces[f]->NN-1; n++)
			{
				int n_ID;
				n_ID	= faces[f]->node[n];

				faces[f]->x[k]	+= nodes[whereis_node[n_ID]]->x[k];
			}

			faces[f]->x[k]	/= faces[f]->NN;

			if (faces[f]->x[k]	!= faces[f]->x[k] || fabs(faces[f]->x[k]) == numeric_limits<double>::infinity())
			{
				Error_message_type	error_message;

				error_message.module_name	="COMMON-GRID";
				error_message.location_primary_name	= "Face";
				error_message.location_primary		= faces[f]->ID;
				error_message.location_secondary_name	= "Direction";
				error_message.location_secondary		= k;
				error_message.message	= " Cannot find center of a face!";
				error_message.print_message();
			}
		}



		// 2. Calculate S (Area/ Volume)
		switch (faces[f]->type)
		{
		case F_LINE:
			faces[f]->S	= length(nodes[whereis_node[faces[f]->node[0]]]->x,
								nodes[whereis_node[faces[f]->node[1]]]->x,
								DIM);
			break;

		case F_TRIANGLE:
			faces[f]->S	= area_triangle(nodes[whereis_node[faces[f]->node[0]]]->x,
										nodes[whereis_node[faces[f]->node[1]]]->x,
										nodes[whereis_node[faces[f]->node[2]]]->x,
										DIM);
			break;

		case F_QUADRILATERAL:
			faces[f]->S	= area_quadrilateral(nodes[whereis_node[faces[f]->node[0]]]->x,
											nodes[whereis_node[faces[f]->node[1]]]->x,
											nodes[whereis_node[faces[f]->node[2]]]->x,
											nodes[whereis_node[faces[f]->node[3]]]->x,
											DIM);
			break;
		}

		if (faces[f]->S < 0.0 || faces[f]->S != faces[f]->S || fabs(faces[f]->S) == numeric_limits<double>::infinity())
		{
			Error_message_type	error_message;
			error_message.module_name	="COMMON-GRID";
			error_message.location_primary_name	= "Face";
			error_message.location_primary		= faces[f]->ID;
			error_message.location_secondary_name	= "NONE";
			error_message.location_secondary		= 0;
			error_message.message	= " Cannot calculate the area of a face!";
			error_message.print_message();
		}




		// 3. Find directional vectors
		int n1, n2, n3, n4;
		double x1, y1, z1;
		double x2, y2, z2;
		double x3, y3, z3;
		double S, length_x1, length_x2;
		double error_ck;

		switch (faces[f]->type)
		{
			case F_LINE:
				n1	= faces[f]->node[0];
				n2	= faces[f]->node[1];

				x1	= nodes[whereis_node[n2]]->x[0]	- nodes[whereis_node[n1]]->x[0];
				y1	= nodes[whereis_node[n2]]->x[1]	- nodes[whereis_node[n1]]->x[1];
				S	= sqrt(x1*x1 + y1*y1);

				faces[f]->n[0][0]	= y1/S;
				faces[f]->n[0][1]	= -x1/S;

				faces[f]->n[1][0]	= -faces[f]->n[0][1];
				faces[f]->n[1][1]	= faces[f]->n[0][0];

				// ERROR CHECK
				error_ck = sqrt(faces[f]->n[0][0]*faces[f]->n[0][0] + faces[f]->n[0][1]*faces[f]->n[0][1]);
				if (error_ck != error_ck || fabs(error_ck) == numeric_limits<double>::infinity() || fabs(error_ck - 1.0) > 1.0e-10)
				{
					Error_message_type	error_message;
					error_message.module_name	="COMMON-GRID";
					error_message.location_primary_name	= "Face";
					error_message.location_primary		= faces[f]->ID;
					error_message.location_secondary_name	= "NONE";
					error_message.location_secondary		= 0;
					error_message.message	= " Cannot find the normal direction vector of a face!";
					error_message.print_message();
				}

				error_ck = sqrt(faces[f]->n[1][0]*faces[f]->n[1][0] + faces[f]->n[1][1]*faces[f]->n[1][1]);
				if (error_ck != error_ck || fabs(error_ck) == numeric_limits<double>::infinity() || fabs(error_ck - 1.0) > 1.0e-10)
				{
					Error_message_type	error_message;
					error_message.module_name	="COMMON-GRID";
					error_message.location_primary_name	= "Face";
					error_message.location_primary		= faces[f]->ID;
					error_message.location_secondary_name	= "NONE";
					error_message.location_secondary		= 0;
					error_message.message	= " Cannot find the tangential direction vector of a face!";
					error_message.print_message();
				}
				break;

			case F_TRIANGLE:
				n1	= faces[f]->node[0];
				n2	= faces[f]->node[1];
				n3	= faces[f]->node[3];

				x1	= nodes[whereis_node[n2]]->x[0]	- nodes[whereis_node[n1]]->x[0];
				y1	= nodes[whereis_node[n2]]->x[1]	- nodes[whereis_node[n1]]->x[1];
				z1	= nodes[whereis_node[n2]]->x[2]	- nodes[whereis_node[n1]]->x[2];
				length_x1	= sqrt(x1*x1 + y1*y1 + z1*z1);

				x2	= nodes[whereis_node[n3]]->x[0]	- nodes[whereis_node[n1]]->x[0];
				y2	= nodes[whereis_node[n3]]->x[1]	- nodes[whereis_node[n1]]->x[1];
				z2	= nodes[whereis_node[n3]]->x[2]	- nodes[whereis_node[n1]]->x[2];

				faces[f]->n[0][0]	=  0.5*(y1*z2 - z1*y2) 	/ faces[f]->S;
				faces[f]->n[0][1]	=  -0.5*(x1*z2 - z1*x2) / faces[f]->S;
				faces[f]->n[0][2]	=  0.5*(x1*y2 - y1*x2) 	/ faces[f]->S;

				faces[f]->n[1][0]	= x1 / length_x1;
				faces[f]->n[1][1]	= y1 / length_x1;
				faces[f]->n[1][2]	= z1 / length_x1;

				faces[f]->n[2][0]	=  faces[f]->n[0][1]*faces[f]->n[1][2] - faces[f]->n[0][2]*faces[f]->n[1][1];
				faces[f]->n[2][1]	= -faces[f]->n[0][0]*faces[f]->n[1][2] + faces[f]->n[0][2]*faces[f]->n[1][0];
				faces[f]->n[2][2]	=  faces[f]->n[0][0]*faces[f]->n[1][1] - faces[f]->n[0][1]*faces[f]->n[1][0];


				// ERROR CHECK
				error_ck = sqrt(faces[f]->n[0][0]*faces[f]->n[0][0] + faces[f]->n[0][1]*faces[f]->n[0][1] + faces[f]->n[0][2]*faces[f]->n[0][2]);
				if (error_ck != error_ck || fabs(error_ck) == numeric_limits<double>::infinity() || fabs(error_ck - 1.0) > 1.0e-10)
				{
					Error_message_type	error_message;
					error_message.module_name	="COMMON-GRID";
					error_message.location_primary_name	= "Face";
					error_message.location_primary		= faces[f]->ID;
					error_message.location_secondary_name	= "NONE";
					error_message.location_secondary		= 0;
					error_message.message	= " Cannot find the normal direction vector of a face!";
					error_message.print_message();
				}

				error_ck = sqrt(faces[f]->n[1][0]*faces[f]->n[1][0] + faces[f]->n[1][1]*faces[f]->n[1][1] + faces[f]->n[1][2]*faces[f]->n[1][2]);
				if (error_ck != error_ck || fabs(error_ck) == numeric_limits<double>::infinity() || fabs(error_ck - 1.0) > 1.0e-10)
				{
					Error_message_type	error_message;
					error_message.module_name	="COMMON-GRID";
					error_message.location_primary_name	= "Face";
					error_message.location_primary		= faces[f]->ID;
					error_message.location_secondary_name	= "NONE";
					error_message.location_secondary		= 0;
					error_message.message	= " Cannot find the tangential direction_1 vector of a face!";
					error_message.print_message();
				}

				error_ck = sqrt(faces[f]->n[2][0]*faces[f]->n[2][0] + faces[f]->n[2][1]*faces[f]->n[2][1] + faces[f]->n[2][2]*faces[f]->n[2][2]);
				if (error_ck != error_ck || fabs(error_ck) == numeric_limits<double>::infinity() || fabs(error_ck - 1.0) > 1.0e-10)
				{
					Error_message_type	error_message;
					error_message.module_name	="COMMON-GRID";
					error_message.location_primary_name	= "Face";
					error_message.location_primary		= faces[f]->ID;
					error_message.location_secondary_name	= "NONE";
					error_message.location_secondary		= 0;
					error_message.message	= " Cannot find the tangential direction_2 vector of a face!";
					error_message.print_message();
				}
				break;

			case F_QUADRILATERAL:
				n1	= faces[f]->node[0];
				n2	= faces[f]->node[1];
				n3	= faces[f]->node[3];
				n4	= faces[f]->node[4];

				x1	= nodes[whereis_node[n2]]->x[0]	- nodes[whereis_node[n1]]->x[0];
				y1	= nodes[whereis_node[n2]]->x[1]	- nodes[whereis_node[n1]]->x[1];
				z1	= nodes[whereis_node[n2]]->x[2]	- nodes[whereis_node[n1]]->x[2];

				x2	= nodes[whereis_node[n3]]->x[0]	- nodes[whereis_node[n1]]->x[0];
				y2	= nodes[whereis_node[n3]]->x[1]	- nodes[whereis_node[n1]]->x[1];
				z2	= nodes[whereis_node[n3]]->x[2]	- nodes[whereis_node[n1]]->x[2];
				length_x2	= sqrt(x2*x2 + y2*y2 + z2*z2);

				x3	= nodes[whereis_node[n4]]->x[0]	- nodes[whereis_node[n1]]->x[0];
				y3	= nodes[whereis_node[n4]]->x[1]	- nodes[whereis_node[n1]]->x[1];
				z3	= nodes[whereis_node[n4]]->x[2]	- nodes[whereis_node[n1]]->x[2];


				faces[f]->n[0][0]	=   0.5*(y1*z2 - z1*y2 + y2*z3 - z2*y3) / faces[f]->S;
				faces[f]->n[0][1]	=  -0.5*(x1*z2 - z1*x2 + x2*z3 - z2*x3) / faces[f]->S;
				faces[f]->n[0][2]	=   0.5*(x1*y2 - y1*x2 + x2*y3 - y2*x3) / faces[f]->S;

				faces[f]->n[1][0]	= x2 / length_x2;
				faces[f]->n[1][1]	= y2 / length_x2;
				faces[f]->n[1][2]	= z2 / length_x2;

				faces[f]->n[2][0]	=  faces[f]->n[0][1]*faces[f]->n[1][2] - faces[f]->n[0][2]*faces[f]->n[1][1];
				faces[f]->n[2][1]	= -faces[f]->n[0][0]*faces[f]->n[1][2] + faces[f]->n[0][2]*faces[f]->n[1][0];
				faces[f]->n[2][2]	=  faces[f]->n[0][0]*faces[f]->n[1][1] - faces[f]->n[0][1]*faces[f]->n[1][0];


				// ERROR CHECK
				error_ck = sqrt(faces[f]->n[0][0]*faces[f]->n[0][0] + faces[f]->n[0][1]*faces[f]->n[0][1] + faces[f]->n[0][2]*faces[f]->n[0][2]);
				if (error_ck != error_ck || fabs(error_ck) == numeric_limits<double>::infinity() || fabs(error_ck - 1.0) > 1.0e-10)
				{
					Error_message_type	error_message;
					error_message.module_name	="COMMON-GRID";
					error_message.location_primary_name	= "Face";
					error_message.location_primary		= faces[f]->ID;
					error_message.location_secondary_name	= "NONE";
					error_message.location_secondary		= 0;
					error_message.message	= " Cannot find the normal direction vector of a face!";
					error_message.print_message();
				}

				error_ck = sqrt(faces[f]->n[1][0]*faces[f]->n[1][0] + faces[f]->n[1][1]*faces[f]->n[1][1] + faces[f]->n[1][2]*faces[f]->n[1][2]);
				if (error_ck != error_ck || fabs(error_ck) == numeric_limits<double>::infinity() || fabs(error_ck - 1.0) > 1.0e-10)
				{
					Error_message_type	error_message;
					error_message.module_name	="COMMON-GRID";
					error_message.location_primary_name	= "Face";
					error_message.location_primary		= faces[f]->ID;
					error_message.location_secondary_name	= "NONE";
					error_message.location_secondary		= 0;
					error_message.message	= " Cannot find the tangential direction_1 vector of a face!";
					error_message.print_message();
				}

				error_ck = sqrt(faces[f]->n[2][0]*faces[f]->n[2][0] + faces[f]->n[2][1]*faces[f]->n[2][1] + faces[f]->n[2][2]*faces[f]->n[2][2]);
				if (error_ck != error_ck || fabs(error_ck) == numeric_limits<double>::infinity() || fabs(error_ck - 1.0) > 1.0e-10)
				{
					Error_message_type	error_message;
					error_message.module_name	="COMMON-GRID";
					error_message.location_primary_name	= "Face";
					error_message.location_primary		= faces[f]->ID;
					error_message.location_secondary_name	= "NONE";
					error_message.location_secondary		= 0;
					error_message.message	= " Cannot find the tangential direction_2 vector of a face!";
					error_message.print_message();
				}
				break;
		}
	}
}

void processing_face_data2(unsigned int DIM, unsigned int NFM, vector<NODE_CLASS>	&nodes, vector<FACE_CLASS>	&faces)
{
	for (int f = 1; f <= NFM; f++)
	{
		// 1. Calculate Xc
		for (int k = 0; k <= DIM-1; k++)
		{
			faces[f].x[k]	= 0.0;

			for (int n	= 0; n <= faces[f].NN-1; n++)	faces[f].x[k]	+= nodes[faces[f].node[n]].x[k];
			faces[f].x[k]	/= faces[f].NN;

			if (faces[f].x[k]	!= faces[f].x[k] || fabs(faces[f].x[k]) == numeric_limits<double>::infinity())
			{
				Error_message_type	error_message;

				error_message.module_name	="COMMON-GRID";
				error_message.location_primary_name	= "Face";
				error_message.location_primary		= faces[f].ID;
				error_message.location_secondary_name	= "Direction";
				error_message.location_secondary		= k;
				error_message.message	= " Cannot find center of a face!";
				error_message.print_message();
			}
		}



		// 2. Calculate S (Area/ Volume)
		switch (faces[f].type)
		{
		case F_LINE:
			faces[f].S	= length(nodes[faces[f].node[0]].x,	nodes[faces[f].node[1]].x, DIM);
			break;

		case F_TRIANGLE:
			faces[f].S	= area_triangle(nodes[faces[f].node[0]].x, nodes[faces[f].node[1]].x, nodes[faces[f].node[2]].x, DIM);
			break;

		case F_QUADRILATERAL:
			faces[f].S	= area_quadrilateral(nodes[faces[f].node[0]].x, nodes[faces[f].node[1]].x, nodes[faces[f].node[2]].x, nodes[faces[f].node[3]].x, DIM);
			break;
		}

		if (faces[f].S < 0.0 || faces[f].S != faces[f].S || fabs(faces[f].S) == numeric_limits<double>::infinity())
		{
			Error_message_type	error_message;
			error_message.module_name	="COMMON-GRID";
			error_message.location_primary_name	= "Face";
			error_message.location_primary		= faces[f].ID;
			error_message.location_secondary_name	= "NONE";
			error_message.location_secondary		= 0;
			error_message.message	= " Cannot calculate the area of a face!";
			error_message.print_message();
		}




		// 3. Find directional vectors
		int n1, n2, n3, n4;
		double x1, y1, z1;
		double x2, y2, z2;
		double x3, y3, z3;
		double S, length_x1, length_x2;
		double error_ck;

		switch (faces[f].type)
		{
			case F_LINE:
				n1	= faces[f].node[0];
				n2	= faces[f].node[1];

				x1	= nodes[n2].x[0]	- nodes[n1].x[0];
				y1	= nodes[n2].x[1]	- nodes[n1].x[1];
				S	= sqrt(x1*x1 + y1*y1);

				faces[f].n[0][0]	= y1/S;
				faces[f].n[0][1]	= -x1/S;

				faces[f].n[1][0]	= -faces[f].n[0][1];
				faces[f].n[1][1]	= faces[f].n[0][0];

				// ERROR CHECK
				error_ck = sqrt(faces[f].n[0][0]*faces[f].n[0][0] + faces[f].n[0][1]*faces[f].n[0][1]);
				if (error_ck != error_ck || fabs(error_ck) == numeric_limits<double>::infinity() || fabs(error_ck - 1.0) > 1.0e-10)
				{
					Error_message_type	error_message;
					error_message.module_name	="COMMON-GRID";
					error_message.location_primary_name	= "Face";
					error_message.location_primary		= faces[f].ID;
					error_message.location_secondary_name	= "NONE";
					error_message.location_secondary		= 0;
					error_message.message	= " Cannot find the normal direction vector of a face!";
					error_message.print_message();
				}

				error_ck = sqrt(faces[f].n[1][0]*faces[f].n[1][0] + faces[f].n[1][1]*faces[f].n[1][1]);
				if (error_ck != error_ck || fabs(error_ck) == numeric_limits<double>::infinity() || fabs(error_ck - 1.0) > 1.0e-10)
				{
					Error_message_type	error_message;
					error_message.module_name	="COMMON-GRID";
					error_message.location_primary_name	= "Face";
					error_message.location_primary		= faces[f].ID;
					error_message.location_secondary_name	= "NONE";
					error_message.location_secondary		= 0;
					error_message.message	= " Cannot find the tangential direction vector of a face!";
					error_message.print_message();
				}
				break;

			case F_TRIANGLE:
				n1	= faces[f].node[0];
				n2	= faces[f].node[1];
				n3	= faces[f].node[3];

				x1	= nodes[n2].x[0]	- nodes[n1].x[0];
				y1	= nodes[n2].x[1]	- nodes[n1].x[1];
				z1	= nodes[n2].x[2]	- nodes[n1].x[2];
				length_x1	= sqrt(x1*x1 + y1*y1 + z1*z1);

				x2	= nodes[n3].x[0]	- nodes[n1].x[0];
				y2	= nodes[n3].x[1]	- nodes[n1].x[1];
				z2	= nodes[n3].x[2]	- nodes[n1].x[2];

				faces[f].n[0][0]	=  0.5*(y1*z2 - z1*y2) 	/ faces[f].S;
				faces[f].n[0][1]	=  -0.5*(x1*z2 - z1*x2) / faces[f].S;
				faces[f].n[0][2]	=  0.5*(x1*y2 - y1*x2) 	/ faces[f].S;

				faces[f].n[1][0]	= x1 / length_x1;
				faces[f].n[1][1]	= y1 / length_x1;
				faces[f].n[1][2]	= z1 / length_x1;

				faces[f].n[2][0]	=  faces[f].n[0][1]*faces[f].n[1][2] - faces[f].n[0][2]*faces[f].n[1][1];
				faces[f].n[2][1]	= -faces[f].n[0][0]*faces[f].n[1][2] + faces[f].n[0][2]*faces[f].n[1][0];
				faces[f].n[2][2]	=  faces[f].n[0][0]*faces[f].n[1][1] - faces[f].n[0][1]*faces[f].n[1][0];


				// ERROR CHECK
				error_ck = sqrt(faces[f].n[0][0]*faces[f].n[0][0] + faces[f].n[0][1]*faces[f].n[0][1] + faces[f].n[0][2]*faces[f].n[0][2]);
				if (error_ck != error_ck || fabs(error_ck) == numeric_limits<double>::infinity() || fabs(error_ck - 1.0) > 1.0e-10)
				{
					Error_message_type	error_message;
					error_message.module_name	="COMMON-GRID";
					error_message.location_primary_name	= "Face";
					error_message.location_primary		= faces[f].ID;
					error_message.location_secondary_name	= "NONE";
					error_message.location_secondary		= 0;
					error_message.message	= " Cannot find the normal direction vector of a face!";
					error_message.print_message();
				}

				error_ck = sqrt(faces[f].n[1][0]*faces[f].n[1][0] + faces[f].n[1][1]*faces[f].n[1][1] + faces[f].n[1][2]*faces[f].n[1][2]);
				if (error_ck != error_ck || fabs(error_ck) == numeric_limits<double>::infinity() || fabs(error_ck - 1.0) > 1.0e-10)
				{
					Error_message_type	error_message;
					error_message.module_name	="COMMON-GRID";
					error_message.location_primary_name	= "Face";
					error_message.location_primary		= faces[f].ID;
					error_message.location_secondary_name	= "NONE";
					error_message.location_secondary		= 0;
					error_message.message	= " Cannot find the tangential direction_1 vector of a face!";
					error_message.print_message();
				}

				error_ck = sqrt(faces[f].n[2][0]*faces[f].n[2][0] + faces[f].n[2][1]*faces[f].n[2][1] + faces[f].n[2][2]*faces[f].n[2][2]);
				if (error_ck != error_ck || fabs(error_ck) == numeric_limits<double>::infinity() || fabs(error_ck - 1.0) > 1.0e-10)
				{
					Error_message_type	error_message;
					error_message.module_name	="COMMON-GRID";
					error_message.location_primary_name	= "Face";
					error_message.location_primary		= faces[f].ID;
					error_message.location_secondary_name	= "NONE";
					error_message.location_secondary		= 0;
					error_message.message	= " Cannot find the tangential direction_2 vector of a face!";
					error_message.print_message();
				}
				break;

			case F_QUADRILATERAL:
				n1	= faces[f].node[0];
				n2	= faces[f].node[1];
				n3	= faces[f].node[3];
				n4	= faces[f].node[4];

				x1	= nodes[n2].x[0]	- nodes[n1].x[0];
				y1	= nodes[n2].x[1]	- nodes[n1].x[1];
				z1	= nodes[n2].x[2]	- nodes[n1].x[2];

				x2	= nodes[n3].x[0]	- nodes[n1].x[0];
				y2	= nodes[n3].x[1]	- nodes[n1].x[1];
				z2	= nodes[n3].x[2]	- nodes[n1].x[2];
				length_x2	= sqrt(x2*x2 + y2*y2 + z2*z2);

				x3	= nodes[n4].x[0]	- nodes[n1].x[0];
				y3	= nodes[n4].x[1]	- nodes[n1].x[1];
				z3	= nodes[n4].x[2]	- nodes[n1].x[2];


				faces[f].n[0][0]	=   0.5*(y1*z2 - z1*y2 + y2*z3 - z2*y3) / faces[f].S;
				faces[f].n[0][1]	=  -0.5*(x1*z2 - z1*x2 + x2*z3 - z2*x3) / faces[f].S;
				faces[f].n[0][2]	=   0.5*(x1*y2 - y1*x2 + x2*y3 - y2*x3) / faces[f].S;

				faces[f].n[1][0]	= x2 / length_x2;
				faces[f].n[1][1]	= y2 / length_x2;
				faces[f].n[1][2]	= z2 / length_x2;

				faces[f].n[2][0]	=  faces[f].n[0][1]*faces[f].n[1][2] - faces[f].n[0][2]*faces[f].n[1][1];
				faces[f].n[2][1]	= -faces[f].n[0][0]*faces[f].n[1][2] + faces[f].n[0][2]*faces[f].n[1][0];
				faces[f].n[2][2]	=  faces[f].n[0][0]*faces[f].n[1][1] - faces[f].n[0][1]*faces[f].n[1][0];


				// ERROR CHECK
				error_ck = sqrt(faces[f].n[0][0]*faces[f].n[0][0] + faces[f].n[0][1]*faces[f].n[0][1] + faces[f].n[0][2]*faces[f].n[0][2]);
				if (error_ck != error_ck || fabs(error_ck) == numeric_limits<double>::infinity() || fabs(error_ck - 1.0) > 1.0e-10)
				{
					Error_message_type	error_message;
					error_message.module_name	="COMMON-GRID";
					error_message.location_primary_name	= "Face";
					error_message.location_primary		= faces[f].ID;
					error_message.location_secondary_name	= "NONE";
					error_message.location_secondary		= 0;
					error_message.message	= " Cannot find the normal direction vector of a face!";
					error_message.print_message();
				}

				error_ck = sqrt(faces[f].n[1][0]*faces[f].n[1][0] + faces[f].n[1][1]*faces[f].n[1][1] + faces[f].n[1][2]*faces[f].n[1][2]);
				if (error_ck != error_ck || fabs(error_ck) == numeric_limits<double>::infinity() || fabs(error_ck - 1.0) > 1.0e-10)
				{
					Error_message_type	error_message;
					error_message.module_name	="COMMON-GRID";
					error_message.location_primary_name	= "Face";
					error_message.location_primary		= faces[f].ID;
					error_message.location_secondary_name	= "NONE";
					error_message.location_secondary		= 0;
					error_message.message	= " Cannot find the tangential direction_1 vector of a face!";
					error_message.print_message();
				}

				error_ck = sqrt(faces[f].n[2][0]*faces[f].n[2][0] + faces[f].n[2][1]*faces[f].n[2][1] + faces[f].n[2][2]*faces[f].n[2][2]);
				if (error_ck != error_ck || fabs(error_ck) == numeric_limits<double>::infinity() || fabs(error_ck - 1.0) > 1.0e-10)
				{
					Error_message_type	error_message;
					error_message.module_name	="COMMON-GRID";
					error_message.location_primary_name	= "Face";
					error_message.location_primary		= faces[f].ID;
					error_message.location_secondary_name	= "NONE";
					error_message.location_secondary		= 0;
					error_message.message	= " Cannot find the tangential direction_2 vector of a face!";
					error_message.print_message();
				}
				break;
		}
	}
}



/*
 * ================================================================================
 * 		Process cell data
 * ================================================================================
 */
void processing_cell_data(unsigned int DIM, unsigned int NFM, unsigned int NCM,
						vector<NODE_CLASS *>	&nodes, vector<int> &whereis_node,
						vector<FACE_CLASS *>	&faces, vector<int> &whereis_face,
						vector<CELL_CLASS *>	&cells, vector<int> &whereis_cell)
{
	// 1. Assign NN and NF
	for (int c	= 0; c <= NCM-1; c++)	cells[c]->assign_NN_NF();



	// 2. Find face and node information of a cell
	int	i, j;
	int	flag;
	int cl, cr;
	vector <int>	filled_nodes(NCM+1, 0);
	vector <int>	face_counter(NCM+1, 0);


	// Step 2.1: FOR ALL FACES, CONSIDERING THE RIGHT CELLS
	for (i = 0; i <= NFM-1; i++)
	{
		int cl	= faces[i]->cl[0];
		int cr	= faces[i]->cr[0];

		int cl_pt	= whereis_cell[cl];
		int cr_pt	= whereis_cell[cr];



		// FOR RIGHT-CELL //
		if (cr >= 1 && cr <= NCM)
		{
			switch (cells[cr_pt]->type)
			{
			// TRIANGLE //
			case C_TRIANGLE:
				if (faces[i]->type == F_LINE)
				{
					if (filled_nodes[cr] == 0)		// FIRST FILLING
					{
						cells[cr_pt]->node[0]	= faces[i]->node[1];
						cells[cr_pt]->node[1]	= faces[i]->node[0];
						filled_nodes[cr] 			= 2;

						cells[cr_pt]->face[0]	= faces[i]->ID;
						face_counter[cr] 		= 1;
					}
					else if(filled_nodes[cr] == 2)	// SECOND FILLEING
					{
						if (faces[i]->node[0] == cells[cr_pt]->node[0])
						{
							cells[cr_pt]->node[2]	= faces[i]->node[1];
							filled_nodes[cr]	= 3;

							cells[cr_pt]->face[2]	= faces[i]->ID;
							face_counter[cr]	= 2;
						}
						else if(faces[i]->node[1] == cells[cr_pt]->node[1])
						{
							cells[cr_pt]->node[2]	= faces[i]->node[0];
							filled_nodes[cr]	= 3;

							cells[cr_pt]->face[1]	= faces[i]->ID;
							face_counter[cr]	= 2;
						}
						else
						{
							Error_message_type	error;
							error.location_primary_name	= "processing_grid.cpp";
							error.message = "PROBLEM DURING CELL MESH CONSTRUCTION";
							error.print_message();
						}
					}
					else if(filled_nodes[cr] == 3)	// LAST FILLING AND CHECK
					{
						if (face_counter[cr] == 2)
						{
							if (cells[cr_pt]->face[1] <= 0)
							{
								cells[cr_pt]->face[1]	= faces[i]->ID;
								face_counter[cr]	= 3;
							}
							else if (cells[cr_pt]->face[2] <= 0)
							{
								cells[cr_pt]->face[2]	= faces[i]->ID;
								face_counter[cr]	= 3;
							}
						}
					}
				}
				else
				{
					Error_message_type	error;
					error.location_primary_name	= "processing_grid.cpp";
					error.message = "PROBLEM DURING CELL MESH CONSTRUCTION! FOR TRIANGLE CELL, ONLY LINE TYPE FACE IS ALLOWEED. PLEASE CHECK THE MESH DATA!!";
					error.print_message();
				}
				break;


			// TETRAHEDRON //
			case C_TETRAHEDRON:
				if (faces[i]->type == F_TRIANGLE)
				{
					if (filled_nodes[cr] == 0)		// FIRST FILLING
					{
						cells[cr_pt]->node[0]	= faces[i]->node[0];
						cells[cr_pt]->node[1]	= faces[i]->node[1];
						cells[cr_pt]->node[2]	= faces[i]->node[2];
						filled_nodes[cr]			= 3;

						cells[cr_pt]->face[0]	= faces[i]->ID;
						face_counter[cr]				= 1;
					}
					else if(filled_nodes[cr] == 3)	// SECOND FILLING
					{
						if (faces[i]->node[0]	!= cells[cr_pt]->node[0] &&
							faces[i]->node[0]	!= cells[cr_pt]->node[1] &&
							faces[i]->node[0]	!= cells[cr_pt]->node[2])
						{
							cells[cr_pt]->node[3]	= faces[i]->node[0];
							filled_nodes[cr]			= 4;

							cells[cr_pt]->face[1]	= faces[i]->ID;
							face_counter[cr]				= 2;
						}
						else if(faces[i]->node[1]	!= cells[cr_pt]->node[0] &&
								faces[i]->node[1]	!= cells[cr_pt]->node[1] &&
								faces[i]->node[1]	!= cells[cr_pt]->node[2])
						{
							cells[cr_pt]->node[3]	= faces[i]->node[1];
							filled_nodes[cr]			= 4;

							cells[cr_pt]->face[1]	= faces[i]->ID;
							face_counter[cr]				= 2;
						}
						else if(faces[i]->node[2]	!= cells[cr_pt]->node[0] &&
								faces[i]->node[2]	!= cells[cr_pt]->node[1] &&
								faces[i]->node[2]	!= cells[cr_pt]->node[2])
						{
							cells[cr_pt]->node[3]	= faces[i]->node[2];
							filled_nodes[cr]			= 4;

							cells[cr_pt]->face[1]	= faces[i]->ID;
							face_counter[cr]		= 2;
						}
					}
					else if(filled_nodes[cr] == 4)	// SECOND FILLING
					{
						if (face_counter[cr] != 4)
						{
							flag = 0;
							for (j = 0; j <= filled_nodes[cr]-1; j++)
							{
								if (cells[cr_pt]->face[j] == faces[i]->ID)	flag = 1;
							}

							if (flag != 1)
							{
								j	= face_counter[cr];
								cells[cr_pt]->face[j] = faces[i]->ID;
								face_counter[cr]	+= 1;
							}
						}
					}
				}
				else
				{
					Error_message_type	error;
					error.location_primary_name	= "processing_grid.cpp";
					error.message = "PROBLEM DURING CELL MESH CONSTRUCTION! FOR TETRAHEDRON CELL, ONLY TRIANGLE TYPE FACE IS ALLOWEED. PLEASE CHECK THE MESH DATA!!";
					error.print_message();
				}
				break;

			case C_QUADRILATERAL:
				// QUADRILATERAL
				if (faces[i]->type == F_LINE)
				{
					if (filled_nodes[cr] == 0)		// FIRST FILLING
					{
						cells[cr_pt]->node[0]	= faces[i]->node[1];
						cells[cr_pt]->node[1]	= faces[i]->node[0];
						filled_nodes[cr]			= 2;

						cells[cr_pt]->face[0]	= faces[i]->ID;
						face_counter[cr]		= 1;
					}
					else if (filled_nodes[cr] == 2)	// SECOND FILLING
					{
						if (faces[i]->node[0] != cells[cr_pt]->node[0] &&
							faces[i]->node[0] != cells[cr_pt]->node[1] &&
							faces[i]->node[1] != cells[cr_pt]->node[0] &&
							faces[i]->node[1] != cells[cr_pt]->node[1])
						{
							cells[cr_pt]->node[2]	= faces[i]->node[1];
							cells[cr_pt]->node[3]	= faces[i]->node[0];
							filled_nodes[cr]		= 4;

							cells[cr_pt]->face[2]	= faces[i]->ID;
							face_counter[cr]		+= 1;
						}
						else if(faces[i]->node[1]	== cells[cr_pt]->node[1] &&
								faces[i]->node[0]	!= cells[cr_pt]->node[2])
						{
							cells[cr_pt]->face[1]	= faces[i]->ID;
							face_counter[cr]		+= 1;
						}
						else if(faces[i]->node[0]	== cells[cr_pt]->node[0] &&
								faces[i]->node[1]	!= cells[cr_pt]->node[3])
						{
							cells[cr_pt]->face[3]	= faces[i]->ID;
							face_counter[cr] 		+= 1;
						}
					}
					else if (filled_nodes[cr] == 4)
					{
						if (face_counter[cr] >= 0 && face_counter[cr] < 4)
						{
							if (faces[i]->node[0]	== cells[cr_pt]->node[2] &&
								faces[i]->node[1]	== cells[cr_pt]->node[1])
							{
								cells[cr_pt]->face[1]	= faces[i]->ID;
								face_counter[cr]		+= 1;
							}
							else if(faces[i]->node[0]	== cells[cr_pt]->node[0] &&
									faces[i]->node[1]	== cells[cr_pt]->node[3])
							{
								cells[cr_pt]->face[3]	= faces[i]->ID;
								face_counter[cr]		+= 1;
							}
						}
					}
				}
				else
				{
					Error_message_type	error;
					error.location_primary_name	= "processing_grid.cpp";
					error.message = "PROBLEM DURING CELL MESH CONSTRUCTION! FOR QUADRILATERAL CELL, ONLY TRIANGLE TYPE FACE IS ALLOWEED. PLEASE CHECK THE MESH DATA!!";
					error.print_message();
				}
				break;

			// HEXAHEDRON
			case C_HEXAHEDRON:
				if (faces[i]->type == F_QUADRILATERAL)
				{
					if (filled_nodes[cr] == 0)
					{
						cells[cr_pt]->node[0]	= faces[i]->node[0];
						cells[cr_pt]->node[1]	= faces[i]->node[1];
						cells[cr_pt]->node[2]	= faces[i]->node[2];
						cells[cr_pt]->node[3]	= faces[i]->node[3];
						filled_nodes[cr]			= 4;

						cells[cr_pt]->face[0]	= faces[i]->ID;
						face_counter[cr]				= 1;
					}
					else if(filled_nodes[cr] < 8)
					{
						if (cells[cr_pt]->node[0]	== faces[i]->node[0])
						{
							if (cells[cr_pt]->node[1]	== faces[i]->node[3])
							{
								cells[cr_pt]->node[4]	= faces[i]->node[1];
								cells[cr_pt]->node[5]	= faces[i]->node[2];
								filled_nodes[cr]			+= 2;
							}
						}
						else if (cells[cr_pt]->node[0]	== faces[i]->node[3])
						{
							if (cells[cr_pt]->node[1]	== faces[i]->node[2])
							{
								cells[cr_pt]->node[4]	= faces[i]->node[0];
								cells[cr_pt]->node[5]	= faces[i]->node[1];
								filled_nodes[cr]			+= 2;
							}
						}
						else if (cells[cr_pt]->node[0]	== faces[i]->node[2])
						{
							if (cells[cr_pt]->node[1]	== faces[i]->node[1])
							{
								cells[cr_pt]->node[4]	= faces[i]->node[3];
								cells[cr_pt]->node[5]	= faces[i]->node[0];
								filled_nodes[cr] 			+= 2;
							}
						}
						else if (cells[cr_pt]->node[0]	== faces[i]->node[1])
						{
							if (cells[cr_pt]->node[1]	== faces[i]->node[0])
							{
								cells[cr_pt]->node[4]	= faces[i]->node[2];
								cells[cr_pt]->node[5]	= faces[i]->node[3];
								filled_nodes[cr]			+= 2;
							}
						}
						else if (cells[cr_pt]->node[3] == faces[i]->node[0])
						{
							if (cells[cr_pt]->node[2]	== faces[i]->node[1])
							{
								cells[cr_pt]->node[6]	= faces[i]->node[2];
								cells[cr_pt]->node[7]	= faces[i]->node[3];
								filled_nodes[cr] 			+= 2;
							}
						}
						else if (cells[cr_pt]->node[3]	== faces[i]->node[3])
						{
							if (cells[cr_pt]->node[2]	== faces[i]->node[0])
							{
								cells[cr_pt]->node[6]	= faces[i]->node[1];
								cells[cr_pt]->node[7]	= faces[i]->node[2];
								filled_nodes[cr]			+= 2;
							}
						}
						else if (cells[cr_pt]->node[3]	== faces[i]->node[2])
						{
							if (cells[cr_pt]->node[2]	== faces[i]->node[3])
							{
								cells[cr_pt]->node[6]	= faces[i]->node[0];
								cells[cr_pt]->node[7]	= faces[i]->node[1];
								filled_nodes[cr]			+= 2;
							}
						}
						else if (cells[cr_pt]->node[3] == faces[i]->node[1])
						{
							if (cells[cr_pt]->node[2]	== faces[i]->node[2])
							{
								cells[cr_pt]->node[6]	= faces[i]->node[3];
								cells[cr_pt]->node[7]	= faces[i]->node[0];
								filled_nodes[cr] 			+= 2;
							}
						}
					}
					else if(filled_nodes[cr] > 8)
					{
						Error_message_type	error;
						error.location_primary_name	= "processing_grid.cpp";
						error.message = "PROBLEM DURING CELL MESH CONSTRUCTION! FOR HEXAHEDRON CELL, NUMBER OF NODE SHOULD BE 8. PLEASE CHECK THE MESH DATA!!";
						error.print_message();
					}
				}
				else
				{
					Error_message_type	error;
					error.location_primary_name	= "processing_grid.cpp";
					error.message = "PROBLEM DURING CELL MESH CONSTRUCTION! FOR HEXAHEDRON CELL, ONLY QUADRILATERAL TYPE FACE IS ALLOWEED. PLEASE CHECK THE MESH DATA!!";
					error.print_message();
				}

				break;

			// PYRAMID
			case C_PYRAMID:
				// Need to add
				break;

			// WEDGE
			case C_WEDGE:
				// Need to add
				break;
			}
		}


		// FOR LEFT-CELL //
		if (cl >= 1 && cl <= NCM)
		{
			switch (cells[cl_pt]->type)
			{
			case C_TRIANGLE:
				if (faces[i]->type == F_LINE)
				{
					if (filled_nodes[cl]	== 0)
					{
						cells[cl_pt]->node[0]	= faces[i]->node[0];
						cells[cl_pt]->node[1]	= faces[i]->node[1];
						filled_nodes[cl] 			= 2;

						cells[cl_pt]->face[0]	= faces[i]->ID;
						face_counter[cl]				= 1;
					}
					else if(filled_nodes[cl]	== 2)
					{
						if (faces[i]->node[0]	== cells[cl_pt]->node[1])
						{
							cells[cl_pt]->node[2]	= faces[i]->node[1];
							filled_nodes[cl] 			= 3;

							cells[cl_pt]->face[1]	= faces[i]->ID;
							face_counter[cl] 				= 2;
						}
						else if(faces[i]->node[1]	== cells[cl_pt]->node[0])
						{
							cells[cl_pt]->node[2]	= faces[i]->node[0];
							filled_nodes[cl]			= 3;

							cells[cl_pt]->face[2] = faces[i]->ID;
							face_counter[cl] = 2;
						}
						else
						{
							Error_message_type	error;
							error.location_primary_name	= "processing_grid.cpp";
							error.message = "PROBLEM DURING CELL MESH CONSTRACTION!";
							error.print_message();
						}
					}
					else if(filled_nodes[cl] == 3)
					{
						if (face_counter[cl] == 2)
						{
							if (cells[cl_pt]->face[1] <= 0)
							{
								cells[cl_pt]->face[1]	= faces[i]->ID;
								face_counter[cl] 			= 3;
							}
							else if (cells[cl_pt]->face[2] <= 0)
							{
								cells[cl_pt]->face[2]	= faces[i]->ID;
								face_counter[cl] 			= 3;
							}
						}
					}
				}
				else
				{
					Error_message_type	error;
					error.location_primary_name	= "processing_grid.cpp";
					error.message = "PROBLEM DURING CELL MESH CONSTRACTION! FOR TRIANGLE CELL, ONLY LINE TYPE FACE IS ALLOWEED. PLEASE CHECK THE MESH DATA!!";
					error.print_message();
				}
				break;
			case C_TETRAHEDRON:
				if (faces[i]->type == F_TRIANGLE)
				{
					if (filled_nodes[cl] == 0)		// FIRST FILLING
					{
						cells[cl_pt]->node[0]	= faces[i]->node[0];
						cells[cl_pt]->node[1]	= faces[i]->node[2];
						cells[cl_pt]->node[2]	= faces[i]->node[1];
						filled_nodes[cl] 			= 3;

						cells[cl_pt]->face[0]	= faces[i]->ID;
						face_counter[cl] 				= 1;
					}
					else if(filled_nodes[cl] == 3)	// SECOND FILLING
					{
						if (faces[i]->node[0]	!= cells[cl_pt]->node[0] &&
							faces[i]->node[0] 	!= cells[cl_pt]->node[1] &&
							faces[i]->node[0] 	!= cells[cl_pt]->node[2])
						{
							cells[cl_pt]->node[3]	= faces[i]->node[0];
							filled_nodes[cl] 			= 4;

							cells[cl_pt]->face[1] 	= faces[i]->ID;
							face_counter[cl] 				= 2;
						}
						else if(faces[i]->node[1] 	!= cells[cl_pt]->node[0] &&
								faces[i]->node[1] 	!= cells[cl_pt]->node[1] &&
								faces[i]->node[1] 	!= cells[cl_pt]->node[2])
						{
							cells[cl_pt]->node[3] 	= faces[i]->node[1];
							filled_nodes[cl] 			= 4;

							cells[cl_pt]->face[1]	= faces[i]->ID;
							face_counter[cl] 				= 2;
						}
						else if(faces[i]->node[2] 	!= cells[cl_pt]->node[0] &&
								faces[i]->node[2] 	!= cells[cl_pt]->node[1] &&
								faces[i]->node[2] 	!= cells[cl_pt]->node[2])
						{
							cells[cl_pt]->node[3] 	= faces[i]->node[2];
							filled_nodes[cl] 			= 4;

							cells[cl_pt]->face[1] 	= faces[i]->ID;
							face_counter[cl]				= 2;
						}
					}
					else if(filled_nodes[cl] == 4)	// SECOND FILLING
					{
						if (face_counter[cl] != 4)
						{
							flag = 0;
							for (j = 0; j <= filled_nodes[cl]-1; j++)
							{
								if (cells[cl_pt]->face[j] == faces[i]->ID) flag = 1;
							}

							if (flag != 1)
							{
								j	= face_counter[cl];
								cells[cl_pt]->face[j]	= faces[i]->ID;
								face_counter[cl] 				+= 1;
							}
						}
					}
				}
				else
				{
					Error_message_type	error;
					error.location_primary_name	= "processing_grid.cpp";
					error.message = "PROBLEM DURING CELL MESH CONSTRACTION! FOR TETRAHEDRON CELL, ONLY TRIANGLE TYPE FACE IS ALLOWEED. PLEASE CHECK THE MESH DATA!!";
					error.print_message();
				}
				break;

			case C_QUADRILATERAL:
				// QUADRILATERAL
				if (faces[i]->type == F_LINE)
				{
					if (filled_nodes[cl] == 0)		// FIRST FILLING
					{
						cells[cl_pt]->node[0]	= faces[i]->node[0];
						cells[cl_pt]->node[1]	= faces[i]->node[1];
						filled_nodes[cl] 			= 2;

						cells[cl_pt]->face[0]	= faces[i]->ID;
						face_counter[cl] 				= 1;
					}
					else if (filled_nodes[cl] == 2)	// SECOND FILLING
					{
						if (faces[i]->node[0] 	!= cells[cl_pt]->node[0] &&
							faces[i]->node[0] 	!= cells[cl_pt]->node[1] &&
							faces[i]->node[1] 	!= cells[cl_pt]->node[0] &&
							faces[i]->node[1] 	!= cells[cl_pt]->node[1])
						{
							cells[cl_pt]->node[2] 	= faces[i]->node[0];
							cells[cl_pt]->node[3] 	= faces[i]->node[1];
							filled_nodes[cl] 			= 4;

							cells[cl_pt]->face[2] 	= faces[i]->ID;
							face_counter[cl] 				+= 1;
						}
						else if(faces[i]->node[0] 	== cells[cl_pt]->node[1] &&
								faces[i]->node[1] 	!= cells[cl_pt]->node[2])
						{
							cells[cl_pt]->face[1] 	= faces[i]->ID;
							face_counter[cl] 				+= 1;
						}
						else if(faces[i]->node[1] 	== cells[cl_pt]->node[0] &&
								faces[i]->node[0] 	!= cells[cl_pt]->node[3])
						{
							cells[cl_pt]->face[3] 	= faces[i]->ID;
							face_counter[cl] 			+= 1;
						}
					}
					else if (filled_nodes[cl] == 4)
					{
						if (face_counter[cl] >= 0 && face_counter[cl] < 4)
						{
							if (faces[i]->node[0] 	== cells[cl_pt]->node[1] &&
								faces[i]->node[1] 	== cells[cl_pt]->node[2])
							{
								cells[cl_pt]->face[1] 	= faces[i]->ID;
								face_counter[cl] 				+= 1;
							}
							else if(faces[i]->node[0] 	== cells[cl_pt]->node[3] &&
									faces[i]->node[1] 	== cells[cl_pt]->node[0])
							{
								cells[cl_pt]->face[3] 	= faces[i]->ID;
								face_counter[cl] 				+= 1;
							}
						}
					}
				}
				else
				{
					Error_message_type	error;
					error.location_primary_name	= "processing_grid.cpp";
					error.message = "PROBLEM DURING CELL MESH CONSTRACTION! FOR QUADRILATERAL CELL, ONLY LINE TYPE FACE IS ALLOWEED. PLEASE CHECK THE MESH DATA!!";
					error.print_message();
				}
				break;

			case C_HEXAHEDRON:
			// HEXAHEDRON
				if (faces[i]->type == F_QUADRILATERAL)
				{
					if (filled_nodes[cl] == 0)		// FIRST FILLING
					{
						cells[cl_pt]->node[0] 	= faces[i]->node[0];
						cells[cl_pt]->node[1] 	= faces[i]->node[3];
						cells[cl_pt]->node[2] 	= faces[i]->node[2];
						cells[cl_pt]->node[3] 	= faces[i]->node[1];
						filled_nodes[cl] 			= 4;

						cells[cr_pt]->face[0] 	= faces[i]->ID;
						face_counter[cr] 				= 1;
					}
					else if(filled_nodes[cl] < 8)
					{
						if (cells[cl_pt]->node[0] == faces[i]->node[0])
						{
							if (cells[cl_pt]->node[1] == faces[i]->node[1])
							{
								cells[cl_pt]->node[4] 	= faces[i]->node[2];
								cells[cl_pt]->node[5] 	= faces[i]->node[3];
								filled_nodes[cl] 			+= 2;
							}
						}
						else if (cells[cl_pt]->node[0] == faces[i]->node[3])
						{
							if (cells[cl_pt]->node[1] == faces[i]->node[0])
							{
								cells[cl_pt]->node[4] 	= faces[i]->node[2];
								cells[cl_pt]->node[5] 	= faces[i]->node[1];
								filled_nodes[cl] 			+= 2;
							}
						}
						else if (cells[cl_pt]->node[0] == faces[i]->node[2])
						{
							if (cells[cl_pt]->node[1] == faces[i]->node[3])
							{
								cells[cl_pt]->node[4] 	= faces[i]->node[1];
								cells[cl_pt]->node[5] 	= faces[i]->node[0];
								filled_nodes[cl] 			+= 2;
							}
						}
						else if (cells[cl_pt]->node[0] == faces[i]->node[1])
						{
							if (cells[cl_pt]->node[1] 	== faces[i]->node[2])
							{
								cells[cl_pt]->node[4] 	= faces[i]->node[0];
								cells[cl_pt]->node[5] 	= faces[i]->node[3];
								filled_nodes[cl] 			+= 2;
							}
						}
						else if (cells[cl_pt]->node[3] == faces[i]->node[0])
						{
							if (cells[cl_pt]->node[2] 	== faces[i]->node[3])
							{
								cells[cl_pt]->node[6] 	= faces[i]->node[2];
								cells[cl_pt]->node[7] 	= faces[i]->node[1];
								filled_nodes[cl] 			+= 2;
							}
						}
						else if (cells[cl_pt]->node[3] == faces[i]->node[3])
						{
							if (cells[cl_pt]->node[2] 	== faces[i]->node[2])
							{
								cells[cl_pt]->node[6]	 = faces[i]->node[1];
								cells[cl_pt]->node[7] 	= faces[i]->node[0];
								filled_nodes[cl] 			+= 2;
							}
						}
						else if (cells[cl_pt]->node[3] == faces[i]->node[2])
						{
							if (cells[cl_pt]->node[2] 	== faces[i]->node[1])
							{
								cells[cl_pt]->node[6] 	= faces[i]->node[3];
								cells[cl_pt]->node[7] 	= faces[i]->node[0];
								filled_nodes[cl] 			+= 2;
							}
						}
						else if (cells[cl_pt]->node[3] == faces[i]->node[1])
						{
							if (cells[cl_pt]->node[2] 	== faces[i]->node[0])
							{
								cells[cl_pt]->node[6] 	= faces[i]->node[3];
								cells[cl_pt]->node[7]	= faces[i]->node[2];
								filled_nodes[cl] 			+= 2;
							}
						}
					}
					else if (filled_nodes[cl] > 8)
					{
						Error_message_type	error;
						error.location_primary_name	= "processing_grid.cpp";
						error.message = "PROBLEM DURING CELL MESH CONSTRACTION! FOR HEXAHEDRON CELL, NUMBER OF NODE SHOULD BE 8. PLEASE CHECK THE MESH DATA!!";
						error.print_message();
					}
				}
				else
				{
					Error_message_type	error;
					error.location_primary_name	= "processing_grid.cpp";
					error.message = "PROBLEM DURING CELL MESH CONSTRACTION! FOR HEXAHEDRON CELL, ONLY QUADRILATERAL TYPE FACE IS ALLOWEED. PLEASE CHECK THE MESH DATA!!";
					error.print_message();
				}
				break;
			}
		} // End left cells
	}


	// 3. Initialize BC in cell
	for (i = 0; i <= NCM-1; i++)
	{
		cells[i]->BC	= BC_INTERIOR;
	}


	// 4. Calculate Cell centers
	for (i = 0; i <= NCM-1; i++)
	{
		for (j = 0; j <= DIM-1; j++)
		{
			cells[i]->x[j]	= 0.0;
			for (int n = 0; n <= cells[i]->NN-1; n++)
			{
				int node_ID	= cells[i]->node[n];
				cells[i]->x[j] += nodes[whereis_node[node_ID]]->x[j];
			}

			cells[i]->x[j] /= cells[i]->NN;

			if (cells[i]->x[j] != cells[i]->x[j] || fabs(cells[i]->x[j]) == numeric_limits<double>::infinity())
			{
				Error_message_type	error_message;
				error_message.module_name	="COMMON-GRID";
				error_message.location_primary_name	= "Cell";
				error_message.location_primary		= cells[i]->ID;
				error_message.location_secondary_name	= "DIMENSTION";
				error_message.location_secondary		= j;
				error_message.message	= " Problem in a cell center calculation!";
				error_message.print_message();
			}
		}
	}


	// 5. Calculate Cell area/volume
	int n1, n2, n3, n4, n5, n6, n7, n8;
	for (i = 0; i <= NCM-1; i++)
	{
		switch(cells[i]->type)
		{
		case C_TRIANGLE:
			n1	= cells[i]->node[0];
			n2	= cells[i]->node[1];
			n3	= cells[i]->node[2];
			cells[i]->S	= area_triangle(nodes[whereis_node[n1]]->x, nodes[whereis_node[n2]]->x, nodes[whereis_node[n3]]->x, DIM);
			break;

		case C_TETRAHEDRON:
			n1	= cells[i]->node[0];
			n2	= cells[i]->node[1];
			n3	= cells[i]->node[2];
			n4	= cells[i]->node[3];

			cells[i]->S	= volume_tetrahedron(nodes[whereis_node[n1]]->x, nodes[whereis_node[n2]]->x, nodes[whereis_node[n3]]->x, nodes[whereis_node[n4]]->x);
			break;

		case C_QUADRILATERAL:
			n1	= cells[i]->node[0];
			n2	= cells[i]->node[1];
			n3	= cells[i]->node[2];
			n4	= cells[i]->node[3];
			cells[i]->S	= area_quadrilateral(nodes[whereis_node[n1]]->x, nodes[whereis_node[n2]]->x, nodes[whereis_node[n3]]->x, nodes[whereis_node[n4]]->x, DIM);
			break;

		case C_HEXAHEDRON:
			n1	= cells[i]->node[0];
			n2	= cells[i]->node[1];
			n3	= cells[i]->node[2];
			n4	= cells[i]->node[3];
			n5	= cells[i]->node[4];
			n6	= cells[i]->node[5];
			n7	= cells[i]->node[6];
			n8	= cells[i]->node[7];

			cells[i]->S	= volume_hexahedron(nodes[whereis_node[n1]]->x,
											nodes[whereis_node[n2]]->x,
											nodes[whereis_node[n3]]->x,
											nodes[whereis_node[n4]]->x,
											nodes[whereis_node[n5]]->x,
											nodes[whereis_node[n6]]->x,
											nodes[whereis_node[n7]]->x,
											nodes[whereis_node[n8]]->x);
			break;

		case C_PYRAMID:
			n1	= cells[i]->node[0];
			n2	= cells[i]->node[1];
			n3	= cells[i]->node[2];
			n4	= cells[i]->node[3];
			n5	= cells[i]->node[4];
			cells[i]->S	=volume_pyramid(nodes[whereis_node[n1]]->x,
										nodes[whereis_node[n2]]->x,
										nodes[whereis_node[n3]]->x,
										nodes[whereis_node[n4]]->x,
										nodes[whereis_node[n5]]->x);
			break;

		case C_WEDGE:
			n1	= cells[i]->node[0];
			n2	= cells[i]->node[1];
			n3	= cells[i]->node[2];
			n4	= cells[i]->node[3];
			n5	= cells[i]->node[4];
			n6	= cells[i]->node[5];
			cells[i]->S	= volume_wedge(nodes[whereis_node[n1]]->x,
										nodes[whereis_node[n2]]->x,
										nodes[whereis_node[n3]]->x,
										nodes[whereis_node[n4]]->x,
										nodes[whereis_node[n5]]->x,
										nodes[whereis_node[n6]]->x);
			break;

		default:
			Error_message_type	error;
			error.location_primary_name	= "processing_grid.cpp";
			error.message = "It is not supported cell-type. Please use a different cell-type";
			error.print_message();
			break;
		}


		if (cells[i]->S < 0.0 || cells[i]->S != cells[i]->S || fabs(cells[i]->S) == numeric_limits<double>::infinity())
		{
			Error_message_type	error_message;
			error_message.module_name	="COMMON-GRID";
			error_message.location_primary_name	= "Cell";
			error_message.location_primary		= cells[i]->ID;
			error_message.location_secondary_name	= "Type";
			error_message.location_secondary		= cells[i]->type;
			error_message.message	= " Problem in cell area calculation!";
			error_message.print_message();
		}
	}



	// 6. Characteristic length calculation
	double length, aux;
	double xf[4][3];

	for (int c = 0; c <= NCM-1; c++)
	{
		// Initialize variable
		for (i = 0; i <= 3; i++)	for (j = 0; j <= 2; j++)	xf[i][j] = 0.0;

		switch(cells[c]->type)
		{
		case C_TRIANGLE:
			length	= cells[c]->S;
			cells[c]->characteristic_length	= sqrt(length);
			break;

		case C_TETRAHEDRON:
			n1	= cells[c]->node[0];
			n2	= cells[c]->node[1];
			n3	= cells[c]->node[2];
			n4	= cells[c]->node[3];

			for (i = 0; i <= DIM-1; i++) xf[0][i] = (nodes[whereis_node[n1]]->x[i] + nodes[whereis_node[n2]]->x[i] + nodes[whereis_node[n3]]->x[i]) / 3.0;
			for (i = 0; i <= DIM-1; i++) xf[1][i] = (nodes[whereis_node[n1]]->x[i] + nodes[whereis_node[n2]]->x[i] + nodes[whereis_node[n4]]->x[i]) / 3.0;
			for (i = 0; i <= DIM-1; i++) xf[2][i] = (nodes[whereis_node[n2]]->x[i] + nodes[whereis_node[n3]]->x[i] + nodes[whereis_node[n4]]->x[i]) / 3.0;
			for (i = 0; i <= DIM-1; i++) xf[3][i] = (nodes[whereis_node[n1]]->x[i] + nodes[whereis_node[n4]]->x[i] + nodes[whereis_node[n3]]->x[i]) / 3.0;

			aux	= 0.0;
			for (i = 0; i <= DIM-1; i++)	aux += pow(cells[c]->x[i] - xf[0][i], 2.0);
			length	= sqrt(aux);

			for (j = 1; j <= 3; j++)
			{
				aux = 0.0;
				for (i = 0; i <= DIM-1; i++)	aux += pow(cells[c]->x[i] - xf[j][i], 2.0);
				aux 	= sqrt(aux);

				length	= fmin(length, aux);
			}
			cells[c]->characteristic_length	= 2.0 *length;
			break;

		case C_QUADRILATERAL:
			n1	= cells[c]->node[0];
			n2	= cells[c]->node[1];
			n3	= cells[c]->node[2];
			n4	= cells[c]->node[3];

			for (i = 0; i <= DIM-1; i++) xf[0][i] = nodes[whereis_node[n2]]->x[i] - nodes[whereis_node[n1]]->x[i];
			for (i = 0; i <= DIM-1; i++) xf[1][i] = nodes[whereis_node[n4]]->x[i] - nodes[whereis_node[n1]]->x[i];

			length	= 0.0;
			aux		= 0.0;
			for (i = 0; i <= DIM-1; i++)
			{
				length	+= pow(xf[0][i], 2.0);
				aux		+= pow(xf[1][i], 2.0);
			}

			length	= fmin(length, aux);
			cells[c]->characteristic_length	= sqrt(length);
			break;

		case C_HEXAHEDRON:
			n1	= cells[c]->node[0];
			n2	= cells[c]->node[1];
			n3	= cells[c]->node[2];
			n4	= cells[c]->node[3];

			for (i = 0; i <= DIM-1; i++) xf[0][i] = nodes[whereis_node[n2]]->x[i] - nodes[whereis_node[n1]]->x[i];
			for (i = 0; i <= DIM-1; i++) xf[1][i] = nodes[whereis_node[n3]]->x[i] - nodes[whereis_node[n1]]->x[i];
			for (i = 0; i <= DIM-1; i++) xf[2][i] = nodes[whereis_node[n4]]->x[i] - nodes[whereis_node[n1]]->x[i];


			length	= 0.0;
			aux		= 0.0;
			for (i = 0; i <= DIM-1; i++)
			{
				length	+= pow(xf[0][i], 2.0);
				aux		+= pow(xf[1][i], 2.0);
			}

			length	= fmin(length, aux);
			aux		= 0.0;
			for (i = 0; i <= DIM-1; i++)
			{
				aux		+= pow(xf[2][i], 2.0);
			}

			length	= fmin(length, aux);
			cells[c]->characteristic_length	= sqrt(length);
			break;

		case C_PYRAMID:
			length	= pow(cells[c]->S, 1.0/3.0);
			cells[c]->characteristic_length	= sqrt(length);
			break;

		case C_WEDGE:
			n1	= cells[c]->node[0];
			n2	= cells[c]->node[1];
			n3	= cells[c]->node[2];
			n4	= cells[c]->node[3];

			for (i = 0; i <= DIM-1; i++) xf[0][i] = nodes[whereis_node[n2]]->x[i] - nodes[whereis_node[n1]]->x[i];
			for (i = 0; i <= DIM-1; i++) xf[1][i] = nodes[whereis_node[n3]]->x[i] - nodes[whereis_node[n1]]->x[i];
			for (i = 0; i <= DIM-1; i++) xf[2][i] = nodes[whereis_node[n4]]->x[i] - nodes[whereis_node[n1]]->x[i];

			aux = 	pow(0.5*(xf[0][1]*xf[1][2] - xf[0][2]*xf[1][1]), 2.0);
			aux	+= 	pow(-0.5*(xf[0][0]*xf[1][2] - xf[0][2]*xf[1][0]), 2.0);
			aux	+=	pow(0.5*(xf[0][0]*xf[1][1] - xf[0][1]*xf[1][0]), 2.0);
			length	= pow(aux, 0.25);

			aux	= 0.0;
			for (i = 0; i <= DIM-1; i++) aux	+= pow(xf[2][i], 2.0);
			aux = sqrt(aux);

			length	= fmin(length, aux);
			cells[c]->characteristic_length	= length;
		break;

		default:
			Error_message_type	error;
			error.location_primary_name	= "processing_grid.cpp";
			error.message = "It is not supported cell-type. Please use a different cell-type";
			error.print_message();
			break;
		}


		if (cells[c]->characteristic_length < 0.0 ||
			cells[c]->characteristic_length != cells[c]->characteristic_length ||
			fabs(cells[c]->characteristic_length) == numeric_limits<double>::infinity())
		{
			Error_message_type	error_message;
			error_message.module_name	="COMMON-GRID";
			error_message.location_primary_name	= "Cell";
			error_message.location_primary		= cells[c]->ID;
			error_message.location_secondary_name	= "Type";
			error_message.location_secondary		= cells[c]->type;
			error_message.message	= " Problem in characteristic length calculation of a cell!";
			error_message.print_message();
		}
	}






}









/*
 * ================================================================================
 * 		Process cell-ghost data
 * ================================================================================
 */
void processing_ghost_cell_data_ver1(unsigned int DIM, unsigned int NFM, unsigned int NCM, unsigned int NGM,
									vector<NODE_CLASS *>	&nodes, vector<int> &whereis_node,
									vector<FACE_CLASS *>	&faces, vector<int> &whereis_face,
									vector<CELL_CLASS *>	&cells, vector<int> &whereis_cell,
									vector<CELL_CLASS *>	&cells_ghost, vector<int> &whereis_cell_ghost)
{
	int counter	= 0;
	int last_pt_cell	= cells.size();
	int ghost_cell_ptr;

	whereis_cell_ghost.resize(NGM+1);

	for (int f = 0; f <= NFM-1; f++)
	{
		if(faces[f]->cr[0]	== 0)	// Boundary face: ghost cell
		{
			counter++;
			faces[f]->cr[0]	= -counter;

			// Create Ghost cell
			ghost_cell_ptr	= (last_pt_cell - 1) + counter;
			cells.push_back(new CELL_CLASS());
			whereis_cell.push_back(0);
			whereis_cell[NCM + counter]	= ghost_cell_ptr;

			cells_ghost.push_back(new CELL_CLASS());


			// 1. Geometry data
			cells[ghost_cell_ptr]->ID	= NCM + counter;
			cells[ghost_cell_ptr]->NN	= faces[f]->NN;
			cells[ghost_cell_ptr]->node.resize(cells[ghost_cell_ptr]->NN);

			for (int j = 0; j <= faces[f]->NN; j++)
			{
				cells[ghost_cell_ptr]->node[j]	= faces[f]->node[j];
			}

			for (int j = 0; j <= DIM-1; j++)
			{
				int cl_ptr	= whereis_cell[faces[f]->cl[0]];
				cells[ghost_cell_ptr]->x[j]	= 2.0*faces[f]->x[j]	- cells[cl_ptr]->x[j];
			}

			cells[ghost_cell_ptr]->S						= faces[f]->S;
			cells[ghost_cell_ptr]->characteristic_length	= faces[f]->S;

			// 2. Connectivity data
			cells[ghost_cell_ptr]->NF		= 1;
			cells[ghost_cell_ptr]->face.resize(cells[ghost_cell_ptr]->NF);

			cells[ghost_cell_ptr]->face[0]	= faces[f]->ID;


			// 3. Grid data
			cells[ghost_cell_ptr]->type	= 0;
			cells[ghost_cell_ptr]->BC	= faces[f]->BC;	// BC

			cells_ghost[counter-1]		= cells[ghost_cell_ptr];
			whereis_cell_ghost[counter]	= counter - 1;
		}
		else if (faces[f]->cl[0] <= 0)	// Need to swap left and right cell
		{
			counter++;
			faces[f]->cl[0]	= faces[f]->cr[0];
			faces[f]->cr[0]	= -counter;

			for (int j = 0; j <= DIM-1; j++)
			{
				for (int k = 0; k <= DIM-1; k++)
				{
					faces[f]->n[j][k]	= -faces[f]->n[j][k];
				}
			}

			// Create Ghost cell
			ghost_cell_ptr	= (last_pt_cell - 1) + counter;
			cells.push_back(new CELL_CLASS());
			whereis_cell.push_back(0);
			whereis_cell[NCM + counter]	= ghost_cell_ptr;

			cells_ghost.push_back(new CELL_CLASS());


			// 1. Geometry data
			cells[ghost_cell_ptr]->ID	= NCM + counter;
			cells[ghost_cell_ptr]->NN	= faces[f]->NN;
			cells[ghost_cell_ptr]->node.resize(cells[ghost_cell_ptr]->NN);

			for (int j = 0; j <= faces[f]->NN; j++)
			{
				cells[ghost_cell_ptr]->node[j]	= faces[f]->node[j];
			}

			for (int j = 0; j <= DIM-1; j++)
			{
				int cl_ptr	= whereis_cell[faces[f]->cl[0]];
				cells[ghost_cell_ptr]->x[j]	= 2.0*faces[f]->x[j]	- cells[cl_ptr]->x[j];
			}

			cells[ghost_cell_ptr]->S						= faces[f]->S;
			cells[ghost_cell_ptr]->characteristic_length	= faces[f]->S;

			// 2. Connectivity data
			cells[ghost_cell_ptr]->NF		= 1;
			cells[ghost_cell_ptr]->face.resize(cells[ghost_cell_ptr]->NF);
			cells[ghost_cell_ptr]->face[0]	= faces[f]->ID;



			// 3. Grid data
			cells[ghost_cell_ptr]->type	= 0;
			cells[ghost_cell_ptr]->BC	= faces[f]->BC;	// BC

			cells_ghost[counter-1]		= cells[ghost_cell_ptr];
			whereis_cell_ghost[counter]	= counter - 1;
		}
	}

	// Error check
	if (counter != NGM)
	{
		Error_message_type	error_message;
		error_message.module_name	="COMMON-GRID: Ghost-cell creation";
		error_message.location_primary_name	= "N/A";
		error_message.location_secondary_name	= "N/A";
		error_message.message	= " Problem in the creation of ghost cells!";
		error_message.print_message();
	}
}





/*
 * ==================================================================
 * 	Calculate distance to wall
 * ==================================================================
 */
void processing_cell_data_distance_to_wall(unsigned int DIM, unsigned int NFM, unsigned int NCM, unsigned int NGM,
											vector<NODE_CLASS *>	&nodes, vector<int> &whereis_node,
											vector<FACE_CLASS *>	&faces, vector<int> &whereis_face,
											vector<CELL_CLASS *>	&cells, vector<int> &whereis_cell,
											vector<CELL_CLASS *>	&cells_ghost, vector<int> &whereis_cell_ghost)
{

	// 1. Initizlize variables
	vector<int>		walk_bound(NCM, 0);
	vector<int>		near_wall(NCM, 0);	// Near wall face ID
	for (int i = 0; i <= NCM-1;	i++) cells[i]->dist_wall	= 1.E10;

	int	n_wall_faces	= 0;
	int	n_walk_cells	= 0;
	int face_ptr;
	int cell_ptr, cl, cr;
	double dist[3];

	// 2. Calculate distance for 1st layer
	for (int i = 0; i <= NGM-1; i++)
	{
		if (cells_ghost[i]->BC	== BC_WALL || cells_ghost[i]->BC	== BC_ANODE || cells_ghost[i]->BC	== BC_CATHODE || cells_ghost[i]->BC	== BC_DIELECTRIC_WALL)
		{
			n_wall_faces++;
			face_ptr	= whereis_face[cells_ghost[i]->face[0]];

			cl			= faces[face_ptr]->cl[0];
			if (cl > 0)
			{
				cell_ptr	= whereis_cell[cl];
			}
			else
			{
				Error_message_type	error_message;
				error_message.module_name	="COMMON-GRID: Calculate-distance_to_wall";
				error_message.location_primary_name		= "FACE ID";
				error_message.location_primary			= cells_ghost[i]->face[0];
				error_message.location_secondary_name	= "N/A";
				error_message.message	= " Problem in the stencil of face!";
				error_message.print_message();
			}

			walk_bound[cell_ptr]	= 1;
			near_wall[cell_ptr]		= cells_ghost[i]->face[0];
			n_walk_cells++;

			for (int k = 0; k <= DIM-1;	k++)	dist[k]	= cells[cell_ptr]->x[k] - faces[face_ptr]->x[k];

            double d = 0.0;
            for (int k = 0; k <= DIM-1;	k++)	d	+= dist[k]*dist[k];
            d	= sqrt(d);

            cells[cell_ptr]->dist_wall	= d;
		}
	}

	// 2. Calculate path
	int n_pass = 0;
	int wall;
	while (n_walk_cells > 0 && n_pass <= MAX_NUM_PATH)
	{
		n_pass++;
        n_walk_cells	= 0;

        for (int i = 0; i <= NCM-1; i++)
    	{

        	if (walk_bound[i] == n_pass)
        	{
        		wall	= near_wall[i];	// Near wall face ID

        		for (int f = 0; f <= cells[i]->NF-1; f++)
        		{
        			face_ptr	= whereis_face[cells[i]->face[f]];
        			cl			= faces[face_ptr]->cl[0]	+ faces[face_ptr]->cr[0]	- cells[i]->ID;

        			if (cl	> 0)
        			{
        				cell_ptr	= whereis_cell[cl];

        				if (walk_bound[cell_ptr]	!= 1)
        				{
        					if (near_wall[cell_ptr]	== 0)
        					{
        						near_wall[cell_ptr]		= wall;
        						walk_bound[cell_ptr]	= n_pass + 1;
        						n_walk_cells++;

        						for (int k = 0; k <= DIM-1;	k++)	dist[k]	= cells[cell_ptr]->x[k] - faces[whereis_face[wall]]->x[k];

        						double d = 0.0;
        						for (int k = 0; k <= DIM-1;	k++)	d	+= dist[k]*dist[k];
        						d =	sqrt(d);

        						cells[cell_ptr]->dist_wall	= d;
        					}
        					else if (near_wall[cell_ptr] != wall)
        					{
        						for (int k = 0; k <= DIM-1;	k++)	dist[k]	= cells[cell_ptr]->x[k] - faces[whereis_face[wall]]->x[k];

								double d = 0.0;
								for (int k = 0; k <= DIM-1;	k++)	d	+= dist[k]*dist[k];
								d	= sqrt(d);

								if (d < cells[cell_ptr]->dist_wall)
								{
									near_wall[cell_ptr]		= wall;
									walk_bound[cell_ptr]	= n_pass+1;
									n_walk_cells++;

									cells[cell_ptr]->dist_wall	= d;
								}
        					}
        				}
        			}
        		}
        	}
    	}
	}


	for (int i = 0; i <= NGM-1;i++)
	{
		face_ptr	= whereis_face[cells_ghost[i]->face[0]];
        cl			= faces[face_ptr]->cl[0];
        cell_ptr	= whereis_cell[cl];

        cells_ghost[i]->dist_wall	= -cells[cell_ptr]->dist_wall;
	}


	for (int f = 0; f <= NFM-1; f++)
	{
		cl = faces[f]->cl[0];
		cr = faces[f]->cr[0];

		if (cr > 0) faces[f]->dist_wall	= 0.5 * (cells[whereis_cell[cl]]->dist_wall + cells[whereis_cell[cr]]->dist_wall);
		else		faces[f]->dist_wall	= 0.5 * (cells[whereis_cell[cl]]->dist_wall + cells[whereis_cell_ghost[-cr]]->dist_wall);
		face_ptr	= whereis_face[near_wall[whereis_cell[cl]]];

		double d = 0.0;
        for (int k = 0; k <= DIM-1; k++)	d += faces[f]->n[0][k]*faces[face_ptr]->n[0][k];
        faces[f]->n_dot_wall	= 1.0 - fabs(d);
	}
}









/*
 * ================================================================================
 * 		Stencil Finder
 * ================================================================================
 */
void processing_stencil_ver1(unsigned int DIM, unsigned int NNM, unsigned int NFM, unsigned int NCM, unsigned int NGM,
							vector<NODE_CLASS *>	&nodes, vector<int> &whereis_node,
							vector<FACE_CLASS *>	&faces, vector<int> &whereis_face,
							vector<CELL_CLASS *>	&cells, vector<int> &whereis_cell,
							vector<CELL_CLASS *>	&cells_ghost, vector<int> &whereis_cell_ghost)
{
	int	i, j, k;

	int cr, cl, cll, crr;
	int cl1, cl2, cr1, cr2;
	int index;
	int max_cell_sharing_node;
	int trial_cell;
	int face_ID;

	CVector_MK xf;					// CELL CENTROID
	CVector_MK xc;					// CELL CENTROID
	CVector_MK xfc;
	CVector_MK normal_vector;		// FACE NORMAL VECTOR
	CVector_MK tan_vector;			// FACE TANGENTIAL VECTOR

	double x, y, z;
	double max_dot, dot;
	double min_dot2, max_dot2, dot2;



	// Step 1: Assign Basic variables
	vector<int> cell_per_node(NNM+1, 0);			// INDICATOR EACH NODES ARE BELONG TO HOW MANY CELLS


	// Step 2: FIND OUT EACH NODES ARE BELONG TO HOW MANY CELLS
	int node_ID;
	for (i = 0; i <= NCM-1; i++)
	{
		for (j = 0; j <= cells[i]->NN-1; j++)
		{
			node_ID	= cells[i]->node[j];
			cell_per_node[node_ID]++;
		}
	}

	for (i = 0; i <= NGM-1; i++)
	{
		for (j = 0; j <= cells_ghost[i]->NN-1; j++)
		{
			node_ID	= cells_ghost[i]->node[j];
			cell_per_node[node_ID]++;
		}
	}

	max_cell_sharing_node	= 0;
	for (i = 1; i <= NNM; i++)
	{
		max_cell_sharing_node = max(cell_per_node[i], max_cell_sharing_node);
	}
	for (i = 0; i <= NNM; i++) cell_per_node[i] = 0;



	// Step 3: GET SHARED CELL ID NUMBER FOR EACH NODES
	vector < vector <int> >cell_sharing_node	= vector_2D(NNM+1, max_cell_sharing_node, 0);	// CELL ID# WHICH IS SHARING THE NODE
	for (i = 0; i <= NCM-1; i++)
	{
		for (j = 0; j <= cells[i]->NN-1; j++)
		{
			node_ID	= cells[i]->node[j];
			index = cell_per_node[node_ID];

			if(node_ID == 0)
			{
				int aa = 0;
				aa = 1;
			}
			cell_sharing_node[node_ID][index] = i;
			cell_per_node[node_ID]++;
		}
	}

	for (i = 0; i <= NNM-1; i++)
	{
		nodes[i]->NSC	= cell_per_node[nodes[i]->ID];
		nodes[i]->C_list.resize(nodes[i]->NSC);
		for (j = 0; j <= nodes[i]->NSC-1; j++)
		{
			nodes[i]->C_list[j]	= cell_sharing_node[nodes[i]->ID][j];
		}
	}



	// Step 4: GET NEIGHBOR CELL ID OF CELL i
	vector <int> neighbor_per_cell (NCM+1, 0);														// NUMBER OF NEIGHBOR OF CELLS
	vector < vector <int> >	cellID_neighbor_cell = vector_2D(NCM+1, 6*max_cell_sharing_node, 0);	// CELL ID NUMBERS OF CELL i

	for (i = 1; i <= NCM; i++)
	{
		vector <int> temp_neighbor(6*max_cell_sharing_node, 0);		// TEMPERARY VALUE FOR FINDING NEIGHBOR CELLS
		neighbor_per_cell[i]	= cells[whereis_cell[i]]->NF;		// NUMBER OF NEIGHBOR CELLS OF CELL i

		for (j = 0; j <= neighbor_per_cell[i]-1; j++)
		{
			face_ID	= cells[whereis_cell[i]]->face[j];

			cr = faces[whereis_face[face_ID]]->cl[0] + faces[whereis_face[face_ID]]->cr[0] - i;	// GET RIGHT CELL ID
			temp_neighbor[j] = cr;
		}

		cellID_neighbor_cell[i].resize(neighbor_per_cell[i]);

		for (j = 0; j <= neighbor_per_cell[i]-1; j++)
		{
			cellID_neighbor_cell[i][j] = temp_neighbor[j];
		}
	}




	// Step 5: Find Stencil
	for (i = 0; i <= NFM-1; i++)
	{
		cl = faces[i]->cl[0];	// LEFT CELL
		cr = faces[i]->cr[0];	// RIGHT CELL



		// Step 5.1: FIND CLL / CL1 / CL2
		cll 	= 0;
		max_dot = -1.0;

		cl1			= 0;
		cl2			= 0;
		max_dot2	= -1.0;
		min_dot2	= 1.0;

		// GET NORMAL/TANGENTIAL VECTOR TOWARD TO CL
		for (k = 0; k <= DIM - 1; k++)
		{
			normal_vector.xyz[k]	= -faces[i]->n[0][k];
			tan_vector.xyz[k]		=  faces[i]->n[0][k];
		}

		// GET CELL CENTROID OF CL
		for (k = 0; k <= DIM - 1; k++)	xf.xyz[k]	= cells[whereis_cell[cl]]->x[k];

		// CALCULATE DISTANCE TO CELL-CENTER OF EACH NEIGHBORING CELL
		for (j = 0; j <= neighbor_per_cell[cl]-1; j++)
		{
			trial_cell = cellID_neighbor_cell[cl][j];

			if (trial_cell > 0)
			{
				for (k = 0; k <= DIM - 1; k++)
					xc.xyz[k]	= cells[whereis_cell[trial_cell]]->x[k];
			}
			else
			{
				for (k = 0; k <= DIM - 1; k++)
					xc.xyz[k]	= cells_ghost[whereis_cell_ghost[-trial_cell]]->x[k];
			}

			xfc	= xc - xf;
			xfc.normalize();

			dot 	= normal_vector % xfc;
			dot2	= tan_vector	% xfc;

			if (dot > max_dot)
			{
				max_dot = dot;
				cll = trial_cell;
			}

			if (dot2 > max_dot2)
			{
				max_dot2	= dot2;
				cl1			= trial_cell;
			}

			if (dot2 < min_dot2)
			{
				min_dot2	= dot2;
				cl2			= trial_cell;
			}
		}

		if (cll == 0)	faces[i]->cl[CLL]	= cl;
		else			faces[i]->cl[CLL]	= cll;

		if (cl1 == 0)	faces[i]->cl[CL1]	= cl;
		else			faces[i]->cl[CL1]	= cl1;

		if (cl2 == 0)	faces[i]->cl[CL2]	= cl;
		else			faces[i]->cl[CL2]	= cl2;



		// Step 5.2: FIND CRR
		if (cr >= 1)	// NOT A GHOST CELL
		{
			crr = 0;
			max_dot = -1.0;

			cr1			= 0;
			cr2			= 0;
			max_dot2	= -1.0;
			min_dot2	= 1.0;

			// GET NORMAL VECTOR TOWARD TO CR
			for (k = 0; k <= DIM - 1; k++)	normal_vector.xyz[k]	= faces[i]->n[0][k];

			// GET CELL CENTROID OF CR
			for (k = 0; k <= DIM-1; k++)	xf.xyz[k]	= cells[whereis_cell[cr]]->x[k];

			// CALCULATE DISTANCE TO CELL-CENTER OF EACH NEIGHBORING CELL
			for (j = 0; j <= neighbor_per_cell[cr]-1; j++)
			{
				trial_cell	= cellID_neighbor_cell[cr][j];

				if (trial_cell > 0)
				{
					for (k = 0; k <= DIM - 1; k++)
						xc.xyz[k]	= cells[whereis_cell[trial_cell]]->x[k];
				}
				else
				{
					for (k = 0; k <= DIM - 1; k++)
						xc.xyz[k]	= cells_ghost[whereis_cell_ghost[-trial_cell]]->x[k];
				}

				xfc	= xc - xf;
				xfc.normalize();

				dot 	= normal_vector % xfc;
				dot2	= tan_vector	% xfc;

				if (dot > max_dot)
				{
					max_dot	= dot;
					crr 	= trial_cell;
				}

				if (dot2 > max_dot2)
				{
					max_dot2	= dot2;
					cr1			= trial_cell;
				}

				if (dot2 < min_dot2)
				{
					min_dot2	= dot2;
					cr2			= trial_cell;
				}
			}

			if (crr == 0)	faces[i]->cr[CRR]	= cr;
			else			faces[i]->cr[CRR]	= crr;

			if (cr1 == 0)	faces[i]->cr[CR1]	= cr;
			else			faces[i]->cr[CR1]	= cr1;

			if (cr2 == 0)	faces[i]->cr[CR2]	= cr;
			else			faces[i]->cr[CR2]	= cr2;
		}
		else
		{
			faces[i]->cr[CRR]	= cr;
			faces[i]->cr[CR1]	= cr;
			faces[i]->cr[CR2]	= cr;
		}
	}
}




void processing_cell_data2(unsigned int DIM, unsigned int NFM, unsigned int NCM, vector<NODE_CLASS>	&nodes, vector<FACE_CLASS>	&faces, vector<CELL_CLASS>	&cells)
{
	// 1. Calculate Cell centers
	for (int i = 1; i <= NCM; i++)
	{
		for (int j = 0; j <= DIM-1; j++)
		{
			cells[i].x[j]	= 0.0;
			for (int n = 0; n <= cells[i].NN-1; n++)	cells[i].x[j] += nodes[cells[i].node[n]].x[j];
			cells[i].x[j] /= cells[i].NN;

			if (cells[i].x[j] != cells[i].x[j] || fabs(cells[i].x[j]) == numeric_limits<double>::infinity())
			{
				Error_message_type	error_message;
				error_message.module_name	="COMMON-GRID";
				error_message.location_primary_name	= "Cell";
				error_message.location_primary		= cells[i].ID;
				error_message.location_secondary_name	= "DIMENSTION";
				error_message.location_secondary		= j;
				error_message.message	= " Problem in a cell center calculation!";
				error_message.print_message();
			}
		}
	}


	// 5. Calculate Cell area/volume
	int n1, n2, n3, n4, n5, n6, n7, n8;
	for (int i = 1; i <= NCM; i++)
	{
		switch(cells[i].type)
		{
		case C_TRIANGLE:
			n1	= cells[i].node[0];
			n2	= cells[i].node[1];
			n3	= cells[i].node[2];
			cells[i].S	= area_triangle(nodes[n1].x, nodes[n2].x, nodes[n3].x, DIM);
			break;

		case C_TETRAHEDRON:
			n1	= cells[i].node[0];
			n2	= cells[i].node[1];
			n3	= cells[i].node[2];
			n4	= cells[i].node[3];

			cells[i].S	= volume_tetrahedron(nodes[n1].x, nodes[n2].x, nodes[n3].x, nodes[n4].x);
			break;

		case C_QUADRILATERAL:
			n1	= cells[i].node[0];
			n2	= cells[i].node[1];
			n3	= cells[i].node[2];
			n4	= cells[i].node[3];
			cells[i].S	= area_quadrilateral(nodes[n1].x, nodes[n2].x, nodes[n3].x, nodes[n4].x, DIM);
			break;

		case C_HEXAHEDRON:
			n1	= cells[i].node[0];
			n2	= cells[i].node[1];
			n3	= cells[i].node[2];
			n4	= cells[i].node[3];
			n5	= cells[i].node[4];
			n6	= cells[i].node[5];
			n7	= cells[i].node[6];
			n8	= cells[i].node[7];
			cells[i].S	= volume_hexahedron(nodes[n1].x, nodes[n2].x, nodes[n3].x, nodes[n4].x, nodes[n5].x, nodes[n6].x, nodes[n7].x, nodes[n8].x);
			break;

		case C_PYRAMID:
			n1	= cells[i].node[0];
			n2	= cells[i].node[1];
			n3	= cells[i].node[2];
			n4	= cells[i].node[3];
			n5	= cells[i].node[4];
			cells[i].S	=volume_pyramid(nodes[n1].x, nodes[n2].x, nodes[n3].x, nodes[n4].x, nodes[n5].x);
			break;

		case C_WEDGE:
			n1	= cells[i].node[0];
			n2	= cells[i].node[1];
			n3	= cells[i].node[2];
			n4	= cells[i].node[3];
			n5	= cells[i].node[4];
			n6	= cells[i].node[5];
			cells[i].S	= volume_wedge(nodes[n1].x, nodes[n2].x, nodes[n3].x, nodes[n4].x, nodes[n5].x, nodes[n6].x);
			break;

		default:
			Error_message_type	error;
			error.location_primary_name	= "processing_grid.cpp";
			error.message = "It is not supported cell-type. Please use a different cell-type";
			error.print_message();
			break;
		}


		if (cells[i].S < 0.0 || cells[i].S != cells[i].S || fabs(cells[i].S) == numeric_limits<double>::infinity())
		{
			Error_message_type	error_message;
			error_message.module_name	="COMMON-GRID";
			error_message.location_primary_name	= "Cell";
			error_message.location_primary		= cells[i].ID;
			error_message.location_secondary_name	= "Type";
			error_message.location_secondary		= cells[i].type;
			error_message.message	= " Problem in cell area calculation!";
			error_message.print_message();
		}
	}
}
