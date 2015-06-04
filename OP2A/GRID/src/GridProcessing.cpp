/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 4, 2015
 *      			Author: Minkwan Kim
 *
 * GridProcessing.cpp
 * 			-  
 *  
 */




#include <vector>
using namespace std;

#include "../include/Grid.hpp"
#include "../include/GridRead.hpp"


namespace OP2A{
namespace GRID{


void Grid::processingNodeData(const double mesh_factor, bool is_axisymmetric)
{
	for (int n = 1; n <= NNM; n++)
	{
		for (int k = 0; k <= ND-1; k++)
		{
			nodes[n].geo.x[k]	/= mesh_factor;
		}
	}

	if (is_axisymmetric == true)
	{
		for (int n = 1; n <= NNM; n++)
		{
			nodes[n].geo.x[1]	+= 1.0E-9;
		}
	}
}



void Grid::processingFaceData()
{
	for (int f = 0; f <= NFM-1; f++)
	{
		// 1. Calculate Xc
		for (int k = 0; k <= ND-1; k++)
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





}
}
