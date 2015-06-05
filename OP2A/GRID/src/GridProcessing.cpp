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

#include <limits>
#include "../include/Grid.hpp"
#include "../include/GridRead.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"
#include "Common/include/Exception_NaNValue.hpp"
#include "Common/include/Exception_NegativeValue.hpp"

#include "Math/include/AreaCalculation.hpp"


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
			faces[f].geo.x[k]	= 0.0;

			for (int n	= 0; n <= faces[f].geo.NN-1; n++)	faces[f].geo.x[k]	+= faces[f].geo.node_list[n]->geo.x[k];
			faces[f].geo.x[k]	/= faces[f].geo.NN;

			if (faces[f].geo.x[k]	!= faces[f].geo.x[k]) throw Common::ExceptionNaNValue (FromHere(), "Nan value for face center location ");
			if (fabs(faces[f].geo.x[k]) == numeric_limits<double>::infinity()) throw Common::ExceptionInfiniteValue (FromHere(), "Inifnite Value for face center location");
		}


		// 2. Calculate S (Area/ Volume)
		switch (faces[f].geo.type)
		{
		case FaceType::f_line:
			faces[f].geo.S	= Math::length(faces[f].geo.node_list[0]->geo.x, faces[f].geo.node_list[1]->geo.x);
			break;

		case FaceType::f_triangle:
			faces[f].geo.S	= Math::CalAreaTriangle(faces[f].geo.node_list[0]->geo.x, faces[f].geo.node_list[1]->geo.x, faces[f].geo.node_list[2]->geo.x);
			break;

		case FaceType::f_quadrilateral:
			faces[f].geo.S	= Math::CalAreaQuadrilateral(faces[f].geo.node_list[0]->geo.x, faces[f].geo.node_list[1]->geo.x, faces[f].geo.node_list[2]->geo.x, faces[f].geo.node_list[3]->geo.x);
			break;
		}



	/*	// 3. Find directional vectors
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
		}*/
	}
}





}
}
