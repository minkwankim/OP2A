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
#include "Common/include/Exception_NoSuchValue.hpp"


#include "Math/include/AreaCalculation.hpp"
#include "Math/include/OP2A_Vector.hpp"



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
	for (int f = 1; f <= NFM; f++)
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



		// 3. Find directional vectors
		Math::VECTOR	face_normal;
		Math::VECTOR 	face_tan;

		if (faces[f].geo.type == FaceType::f_line)
		{
			Math::VECTOR	n12(faces[f].geo.node_list[0]->geo.x, faces[f].geo.node_list[1]->geo.x);
			n12.normalize();

			face_normal = n12.rotate(-MATH_PI/2, Math::VectorDirection::VectorDirection_Z);


			faces[f].geo.n[0][0]	= face_normal(1);
			faces[f].geo.n[0][1]	= face_normal(2);

			faces[f].geo.n[1][0]	= n12(1);
			faces[f].geo.n[1][1]	= n12(2);
		}
		else if(faces[f].geo.type == FaceType::f_triangle)
		{
			Math::VECTOR X1(faces[f].geo.node_list[1]->geo.x, faces[f].geo.node_list[0]->geo.x);
			Math::VECTOR X2(faces[f].geo.node_list[2]->geo.x, faces[f].geo.node_list[0]->geo.x);

			face_normal	= VectorCrossProduct(X1, X2);

			face_normal.normalize();
			X1.normalize();

			face_tan	= VectorCrossProduct(face_normal, X1);
			face_tan.normalize();

			faces[f].geo.n[0][0]	=  face_normal(1);
			faces[f].geo.n[0][1]	=  face_normal(2);
			faces[f].geo.n[0][2]	=  face_normal(3);

			faces[f].geo.n[1][0]	=  X1(1);
			faces[f].geo.n[1][1]	=  X1(2);
			faces[f].geo.n[1][2]	=  X1(3);

			faces[f].geo.n[2][0]	=  face_tan(1);
			faces[f].geo.n[2][1]	=  face_tan(2);
			faces[f].geo.n[2][2]	=  face_tan(3);
		}
		else if(faces[f].geo.type == FaceType::f_quadrilateral)
		{
			Math::VECTOR V1(faces[f].geo.node_list[1]->geo.x, faces[f].geo.node_list[0]->geo.x);
			Math::VECTOR V2(faces[f].geo.node_list[2]->geo.x, faces[f].geo.node_list[0]->geo.x);
			Math::VECTOR V3(faces[f].geo.node_list[3]->geo.x, faces[f].geo.node_list[0]->geo.x);

			face_normal	= NormalFromThreePoint(V1, V2, V3);
			face_normal.normalize();

			V2.normalize();

			face_tan	= VectorCrossProduct(face_normal, V2);
			face_tan.normalize();

			faces[f].geo.n[0][0]	=  face_normal(1);
			faces[f].geo.n[0][1]	=  face_normal(2);
			faces[f].geo.n[0][2]	=  face_normal(3);

			faces[f].geo.n[1][0]	=  V2(1);
			faces[f].geo.n[1][1]	=  V2(2);
			faces[f].geo.n[1][2]	=  V2(3);

			faces[f].geo.n[2][0]	=  face_tan(1);
			faces[f].geo.n[2][1]	=  face_tan(2);
			faces[f].geo.n[2][2]	=  face_tan(3);
		}
		else
		{
			throw Common::ExceptionNoSuchValue (FromHere(), "Selected Face type is not supported");
		}

	}
}





}
}
