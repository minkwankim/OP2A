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
#include "../include/Exception_CellConstruction.hpp"
#include "../include/Exception_GridDataMismatch.hpp"


#include "Common/include/Exception_InfiniteValue.hpp"
#include "Common/include/Exception_NaNValue.hpp"
#include "Common/include/Exception_NegativeValue.hpp"
#include "Common/include/Exception_NoSuchValue.hpp"


#include "Math/include/AreaCalculation.hpp"
#include "Math/include/MathMisc.hpp"
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
			nodes[n].geo.x[1]	+= 1.0E-8;
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
			if (fabs(faces[f].geo.x[k]) == numeric_limits<double>::infinity()) throw Common::ExceptionInfiniteValue (FromHere(), "Infinite Value for face center location");
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

			//face_normal(1) = n12.rotate(-MATH_PI/2, Math::VectorDirection::VectorDirection_Z);
			face_normal(1) = n12(2);
			face_normal(2) = -n12(1);



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


void Grid::processingCellData()
{
	// 2. Find face and node information of a cell
	vector <int>	filled_nodes(NCM+1, 0);
	vector <int>	face_counter(NCM+1, 0);

	for (int f = 1; f <= NFM; f++)
	{
		if (faces[f].geo.cr[0] != NULL)
		{
			switch(faces[f].geo.cr[0]->geo.type)
			{
			case CellType::triangle:
				if (faces[f].geo.type == FaceType::f_line)
				{
					if (filled_nodes[faces[f].geo.cr[0]->geo.ID] == 0)		// FIRST FILLING
					{
						faces[f].geo.cr[0]->geo.node_list[0]	= faces[f].geo.node_list[1];
						faces[f].geo.cr[0]->geo.node_list[1]	= faces[f].geo.node_list[0];
						filled_nodes[faces[f].geo.cr[0]->geo.ID] = 2;

						faces[f].geo.cr[0]->geo.face_list[0]	= &faces[f];
						face_counter[faces[f].geo.cr[0]->geo.ID] = 1;
					}
					else if(filled_nodes[faces[f].geo.cr[0]->geo.ID] == 2)	// SECOND FILLEING
					{
						if (faces[f].geo.node_list[0] == faces[f].geo.cr[0]->geo.node_list[0])
						{
							faces[f].geo.cr[0]->geo.node_list[2]	= faces[f].geo.node_list[1];
							filled_nodes[faces[f].geo.cr[0]->geo.ID] = 3;

							faces[f].geo.cr[0]->geo.face_list[2]	= &faces[f];
							face_counter[faces[f].geo.cr[0]->geo.ID] = 2;
						}
						else if(faces[f].geo.node_list[1] == faces[f].geo.cr[0]->geo.node_list[1])
						{
							faces[f].geo.cr[0]->geo.node_list[2]	= faces[f].geo.node_list[0];
							filled_nodes[faces[f].geo.cr[0]->geo.ID] = 3;

							faces[f].geo.cr[0]->geo.face_list[1]	= &faces[f];
							face_counter[faces[f].geo.cr[0]->geo.ID] = 2;
						}
						else
						{
							throw ExceptionCellConstruction (FromHere(), "PROBLEM DURING CELL MESH CONSTRUCTION: Triangular cell");
						}
					}
					else if(filled_nodes[faces[f].geo.cr[0]->geo.ID] == 3)	// LAST FILLING AND CHECK
					{
						if (face_counter[faces[f].geo.cr[0]->geo.ID] == 2)
						{
							if (faces[f].geo.cr[0]->geo.face_list[1] == NULL)
							{
								faces[f].geo.cr[0]->geo.face_list[1] = &faces[f];
								face_counter[faces[f].geo.cr[0]->geo.ID] =3;
							}
							else if (faces[f].geo.cr[0]->geo.face_list[2] == NULL)
							{
								faces[f].geo.cr[0]->geo.face_list[2] = &faces[f];
								face_counter[faces[f].geo.cr[0]->geo.ID] =3;
							}
						}
					}
				}
				else
				{
					throw ExceptionCellConstruction (FromHere(), "PROBLEM DURING CELL MESH CONSTRUCTION: Triangular cell --ONLY LINE TYPE FACE IS ALLOWEED. PLEASE CHECK THE MESH DATA");
				}
				break;

			case CellType::tetrahedron:
				if (faces[f].geo.type == FaceType::f_triangle)
				{
					if (filled_nodes[faces[f].geo.cr[0]->geo.ID] == 0)		// FIRST FILLING
					{
						faces[f].geo.cr[0]->geo.node_list[0]	= faces[f].geo.node_list[0];
						faces[f].geo.cr[0]->geo.node_list[1]	= faces[f].geo.node_list[1];
						faces[f].geo.cr[0]->geo.node_list[2]	= faces[f].geo.node_list[2];
						filled_nodes[faces[f].geo.cr[0]->geo.ID] = 3;

						faces[f].geo.cr[0]->geo.face_list[0]	= &faces[f];
						face_counter[faces[f].geo.cr[0]->geo.ID] =1;
					}
					else if(filled_nodes[faces[f].geo.cr[0]->geo.ID] == 3)	// SECOND FILLING
					{
						if (faces[f].geo.node_list[0]	!= faces[f].geo.cr[0]->geo.node_list[0] &&
							faces[f].geo.node_list[0]	!= faces[f].geo.cr[0]->geo.node_list[1] &&
							faces[f].geo.node_list[0]	!= faces[f].geo.cr[0]->geo.node_list[2])
						{
							faces[f].geo.cr[0]->geo.node_list[3]	= faces[f].geo.node_list[0];
							filled_nodes[faces[f].geo.cr[0]->geo.ID] = 4;

							faces[f].geo.cr[0]->geo.face_list[1]	= &faces[f];
							face_counter[faces[f].geo.cr[0]->geo.ID] =2;
						}
						else if(faces[f].geo.node_list[1]	!= faces[f].geo.cr[0]->geo.node_list[0] &&
								faces[f].geo.node_list[1]	!= faces[f].geo.cr[0]->geo.node_list[1] &&
								faces[f].geo.node_list[1]	!= faces[f].geo.cr[0]->geo.node_list[2])
						{
							faces[f].geo.cr[0]->geo.node_list[3]	= faces[f].geo.node_list[1];
							filled_nodes[faces[f].geo.cr[0]->geo.ID] = 4;

							faces[f].geo.cr[0]->geo.face_list[1]	= &faces[f];
							face_counter[faces[f].geo.cr[0]->geo.ID] =2;
						}
						else if(faces[f].geo.node_list[2]	!= faces[f].geo.cr[0]->geo.node_list[0] &&
								faces[f].geo.node_list[2]	!= faces[f].geo.cr[0]->geo.node_list[1] &&
								faces[f].geo.node_list[2]	!= faces[f].geo.cr[0]->geo.node_list[2])
						{
							faces[f].geo.cr[0]->geo.node_list[3]	= faces[f].geo.node_list[2];
							filled_nodes[faces[f].geo.cr[0]->geo.ID] = 4;

							faces[f].geo.cr[0]->geo.face_list[1]	= &faces[f];
							face_counter[faces[f].geo.cr[0]->geo.ID] = 2;
						}
					}
					else if(filled_nodes[faces[f].geo.cr[0]->geo.ID] == 4)	// SECOND FILLING
					{
						if (face_counter[faces[f].geo.cr[0]->geo.ID] != 4)
						{
							bool flag = false;
							for (int j = 0; j <= face_counter[faces[f].geo.cr[0]->geo.ID]-1; j++)
							{
								if (faces[f].geo.cr[0]->geo.face_list[j] == &faces[j])	flag = true;
							}

							if (flag != true)
							{
								int j	= face_counter[faces[f].geo.cr[0]->geo.ID];
								faces[f].geo.cr[0]->geo.face_list[j] = &faces[j];
								face_counter[faces[f].geo.cr[0]->geo.ID]	+= 1;
							}
						}
					}
				}
				else
				{
					throw ExceptionCellConstruction (FromHere(), "PROBLEM DURING CELL MESH CONSTRUCTION: Tetrajedrpm cell --ONLY TRIANGLE TYPE FACE IS ALLOWEED. PLEASE CHECK THE MESH DATA");
				}
				break;

			case CellType::quadrilateral:
				// QUADRILATERAL
				if (faces[f].geo.type == FaceType::f_line)
				{
					if (filled_nodes[faces[f].geo.cr[0]->geo.ID] == 0)		// FIRST FILLING
					{
						faces[f].geo.cr[0]->geo.node_list[0]	= faces[f].geo.node_list[1];
						faces[f].geo.cr[0]->geo.node_list[1]	= faces[f].geo.node_list[0];
						filled_nodes[faces[f].geo.cr[0]->geo.ID] = 2;

						faces[f].geo.cr[0]->geo.face_list[0]	= &faces[f];
						face_counter[faces[f].geo.cr[0]->geo.ID] = 1;
					}
					else if (filled_nodes[faces[f].geo.cr[0]->geo.ID] == 2)	// SECOND FILLING
					{
						if (faces[f].geo.node_list[0] != faces[f].geo.cr[0]->geo.node_list[0] &&
							faces[f].geo.node_list[0] != faces[f].geo.cr[0]->geo.node_list[1] &&
							faces[f].geo.node_list[1] != faces[f].geo.cr[0]->geo.node_list[0] &&
							faces[f].geo.node_list[1] != faces[f].geo.cr[0]->geo.node_list[1])
						{
							faces[f].geo.cr[0]->geo.node_list[2]	= faces[f].geo.node_list[1];
							faces[f].geo.cr[0]->geo.node_list[3]	= faces[f].geo.node_list[0];
							filled_nodes[faces[f].geo.cr[0]->geo.ID] = 4;

							faces[f].geo.cr[0]->geo.face_list[2]	= &faces[f];
							face_counter[faces[f].geo.cr[0]->geo.ID] += 1;
						}
						else if(faces[f].geo.node_list[1] == faces[f].geo.cr[0]->geo.node_list[1] &&
								faces[f].geo.node_list[0] != faces[f].geo.cr[0]->geo.node_list[2])
						{
							faces[f].geo.cr[0]->geo.face_list[1]	= &faces[f];
							face_counter[faces[f].geo.cr[0]->geo.ID] += 1;
						}
						else if(faces[f].geo.node_list[0] == faces[f].geo.cr[0]->geo.node_list[0] &&
								faces[f].geo.node_list[1] != faces[f].geo.cr[0]->geo.node_list[3])
						{
							faces[f].geo.cr[0]->geo.face_list[3]	= &faces[f];
							face_counter[faces[f].geo.cr[0]->geo.ID] += 1;
						}
					}
					else if (filled_nodes[faces[f].geo.cr[0]->geo.ID] == 4)
					{
						if (face_counter[faces[f].geo.cr[0]->geo.ID] >= 0 && face_counter[faces[f].geo.cr[0]->geo.ID] < 4)
						{
							if (faces[f].geo.node_list[0]	== faces[f].geo.cr[0]->geo.node_list[2] &&
								faces[f].geo.node_list[1]	== faces[f].geo.cr[0]->geo.node_list[1])
							{
								faces[f].geo.cr[0]->geo.face_list[1]	= &faces[f];
								face_counter[faces[f].geo.cr[0]->geo.ID] += 1;
							}
							else if(faces[f].geo.node_list[0]	== faces[f].geo.cr[0]->geo.node_list[0] &&
									faces[f].geo.node_list[1]	== faces[f].geo.cr[0]->geo.node_list[3])
							{
								faces[f].geo.cr[0]->geo.face_list[3]	= &faces[f];
								face_counter[faces[f].geo.cr[0]->geo.ID] += 1;
							}
						}
					}
				}
				else
				{
					throw ExceptionCellConstruction (FromHere(), "PROBLEM DURING CELL MESH CONSTRUCTION: QUADRILATERAL cell --ONLY Line TYPE FACE IS ALLOWEED. PLEASE CHECK THE MESH DATA");
				}
				break;

			case CellType::hexahedron:
				if (faces[f].geo.type == FaceType::f_quadrilateral)
				{
					if (filled_nodes[faces[f].geo.cr[0]->geo.ID] == 0)
					{
						faces[f].geo.cr[0]->geo.node_list[0]	= faces[f].geo.node_list[0];
						faces[f].geo.cr[0]->geo.node_list[1]	= faces[f].geo.node_list[1];
						faces[f].geo.cr[0]->geo.node_list[2]	= faces[f].geo.node_list[2];
						faces[f].geo.cr[0]->geo.node_list[3]	= faces[f].geo.node_list[3];
						filled_nodes[faces[f].geo.cr[0]->geo.ID] = 4;

						faces[f].geo.cr[0]->geo.face_list[0]	= &faces[f];
						face_counter[faces[f].geo.cr[0]->geo.ID] = 1;
					}
					else if(filled_nodes[faces[f].geo.cr[0]->geo.ID] < 8)
					{
						if (faces[f].geo.cr[0]->geo.node_list[0]	== faces[f].geo.node_list[0])
						{
							if (faces[f].geo.cr[0]->geo.node_list[1]	== faces[f].geo.node_list[3])
							{
								faces[f].geo.cr[0]->geo.node_list[4]	= faces[f].geo.node_list[1];
								faces[f].geo.cr[0]->geo.node_list[5]	= faces[f].geo.node_list[2];
								filled_nodes[faces[f].geo.cr[0]->geo.ID]			+= 2;
							}
						}
						else if (faces[f].geo.cr[0]->geo.node_list[0]	== faces[f].geo.node_list[3])
						{
							if (faces[f].geo.cr[0]->geo.node_list[1]	== faces[f].geo.node_list[2])
							{
								faces[f].geo.cr[0]->geo.node_list[4]	= faces[f].geo.node_list[0];
								faces[f].geo.cr[0]->geo.node_list[5]	= faces[f].geo.node_list[1];
								filled_nodes[faces[f].geo.cr[0]->geo.ID]			+= 2;
							}
						}
						else if (faces[f].geo.cr[0]->geo.node_list[0]	== faces[f].geo.node_list[2])
						{
							if (faces[f].geo.cr[0]->geo.node_list[1]	== faces[f].geo.node_list[1])
							{
								faces[f].geo.cr[0]->geo.node_list[4]	= faces[f].geo.node_list[3];
								faces[f].geo.cr[0]->geo.node_list[5]	= faces[f].geo.node_list[0];
								filled_nodes[faces[f].geo.cr[0]->geo.ID] 			+= 2;
							}
						}
						else if (faces[f].geo.cr[0]->geo.node_list[0]	== faces[f].geo.node_list[1])
						{
							if (faces[f].geo.cr[0]->geo.node_list[1]	== faces[f].geo.node_list[0])
							{
								faces[f].geo.cr[0]->geo.node_list[4]	= faces[f].geo.node_list[2];
								faces[f].geo.cr[0]->geo.node_list[5]	= faces[f].geo.node_list[3];
								filled_nodes[faces[f].geo.cr[0]->geo.ID]			+= 2;
							}
						}
						else if (faces[f].geo.cr[0]->geo.node_list[3] == faces[f].geo.node_list[0])
						{
							if (faces[f].geo.cr[0]->geo.node_list[2]	== faces[f].geo.node_list[1])
							{
								faces[f].geo.cr[0]->geo.node_list[6]	= faces[f].geo.node_list[2];
								faces[f].geo.cr[0]->geo.node_list[7]	= faces[f].geo.node_list[3];
								filled_nodes[faces[f].geo.cr[0]->geo.ID] 			+= 2;
							}
						}
						else if (faces[f].geo.cr[0]->geo.node_list[3]	== faces[f].geo.node_list[3])
						{
							if (faces[f].geo.cr[0]->geo.node_list[2]	== faces[f].geo.node_list[0])
							{
								faces[f].geo.cr[0]->geo.node_list[6]	= faces[f].geo.node_list[1];
								faces[f].geo.cr[0]->geo.node_list[7]	= faces[f].geo.node_list[2];
								filled_nodes[faces[f].geo.cr[0]->geo.ID]			+= 2;
							}
						}
						else if (faces[f].geo.cr[0]->geo.node_list[3]	== faces[f].geo.node_list[2])
						{
							if (faces[f].geo.cr[0]->geo.node_list[2]	== faces[f].geo.node_list[3])
							{
								faces[f].geo.cr[0]->geo.node_list[6]	= faces[f].geo.node_list[0];
								faces[f].geo.cr[0]->geo.node_list[7]	= faces[f].geo.node_list[1];
								filled_nodes[faces[f].geo.cr[0]->geo.ID]			+= 2;
							}
						}
						else if (faces[f].geo.cr[0]->geo.node_list[3] == faces[f].geo.node_list[1])
						{
							if (faces[f].geo.cr[0]->geo.node_list[2]	== faces[f].geo.node_list[2])
							{
								faces[f].geo.cr[0]->geo.node_list[6]	= faces[f].geo.node_list[3];
								faces[f].geo.cr[0]->geo.node_list[7]	= faces[f].geo.node_list[0];
								filled_nodes[faces[f].geo.cr[0]->geo.ID] 			+= 2;
							}
						}
					}
					else if(filled_nodes[faces[f].geo.cr[0]->geo.ID] > 8)
					{
						throw ExceptionCellConstruction (FromHere(), "PROBLEM DURING CELL MESH CONSTRUCTION! FOR HEXAHEDRON CELL, NUMBER OF NODE SHOULD BE 8. PLEASE CHECK THE MESH DATA!!");
					}
				}
				else
				{
					throw ExceptionCellConstruction (FromHere(), "PROBLEM DURING CELL MESH CONSTRUCTION! FOR HEXAHEDRON CELL, ONLY QUADRILATERAL TYPE FACE IS ALLOWEED. PLEASE CHECK THE MESH DATA!!");
				}
				break;

			case CellType::pyramid:
				break;

			case CellType::wedge:
				break;

			}
		}


		// For left-cells
		if (faces[f].geo.cl[0] != NULL)
		{
			switch(faces[f].geo.cl[0]->geo.type)
			{
			case CellType::triangle:
				if (faces[f].geo.type == FaceType::f_line)
				{
					if (filled_nodes[faces[f].geo.cl[0]->geo.ID] == 0)		// FIRST FILLING
					{
						faces[f].geo.cl[0]->geo.node_list[0]	= faces[f].geo.node_list[0];
						faces[f].geo.cl[0]->geo.node_list[1]	= faces[f].geo.node_list[1];
						filled_nodes[faces[f].geo.cl[0]->geo.ID] = 2;

						faces[f].geo.cl[0]->geo.face_list[0]	= &faces[f];
						face_counter[faces[f].geo.cl[0]->geo.ID] = 1;
					}
					else if(filled_nodes[faces[f].geo.cl[0]->geo.ID] == 2)	// SECOND FILLING
					{
						if (faces[f].geo.node_list[0] == faces[f].geo.cl[0]->geo.node_list[1])
						{
							faces[f].geo.cl[0]->geo.node_list[2]	= faces[f].geo.node_list[1];
							filled_nodes[faces[f].geo.cl[0]->geo.ID] = 3;

							faces[f].geo.cl[0]->geo.face_list[1]	= &faces[f];
							face_counter[faces[f].geo.cl[0]->geo.ID] = 2;
						}
						else if(faces[f].geo.node_list[1] == faces[f].geo.cl[0]->geo.node_list[0])
						{
							faces[f].geo.cl[0]->geo.node_list[2]	= faces[f].geo.node_list[0];
							filled_nodes[faces[f].geo.cl[0]->geo.ID] = 3;

							faces[f].geo.cl[0]->geo.face_list[2]	= &faces[f];
							face_counter[faces[f].geo.cl[0]->geo.ID] = 2;
						}
						else
						{
							throw ExceptionCellConstruction (FromHere(), "PROBLEM DURING CELL MESH CONSTRUCTION: Triangular cell");
						}
					}
					else if(filled_nodes[faces[f].geo.cl[0]->geo.ID] == 3)	// LAST FILLING AND CHECK
					{
						if (face_counter[faces[f].geo.cl[0]->geo.ID] == 2)
						{
							if (faces[f].geo.cl[0]->geo.face_list[1] == NULL)
							{
								faces[f].geo.cl[0]->geo.face_list[1] = &faces[f];
								face_counter[faces[f].geo.cl[0]->geo.ID] =3;
							}
							else if (faces[f].geo.cl[0]->geo.face_list[2] == NULL)
							{
								faces[f].geo.cl[0]->geo.face_list[2] = &faces[f];
								face_counter[faces[f].geo.cl[0]->geo.ID] =3;
							}
						}
					}
				}
				else
				{
					throw ExceptionCellConstruction (FromHere(), "PROBLEM DURING CELL MESH CONSTRUCTION: Triangular cell --ONLY LINE TYPE FACE IS ALLOWEED. PLEASE CHECK THE MESH DATA");
				}
				break;

			case CellType::tetrahedron:
				if (faces[f].geo.type == FaceType::f_triangle)
				{
					if (filled_nodes[faces[f].geo.cl[0]->geo.ID] == 0)		// FIRST FILLING
					{
						faces[f].geo.cl[0]->geo.node_list[0]	= faces[f].geo.node_list[0];
						faces[f].geo.cl[0]->geo.node_list[1]	= faces[f].geo.node_list[2];
						faces[f].geo.cl[0]->geo.node_list[2]	= faces[f].geo.node_list[1];
						filled_nodes[faces[f].geo.cl[0]->geo.ID] = 3;

						faces[f].geo.cl[0]->geo.face_list[0]	= &faces[f];
						face_counter[faces[f].geo.cl[0]->geo.ID] =1;
					}
					else if(filled_nodes[faces[f].geo.cl[0]->geo.ID] == 3)	// SECOND FILLING
					{
						if (faces[f].geo.node_list[0]	!= faces[f].geo.cl[0]->geo.node_list[0] &&
							faces[f].geo.node_list[0]	!= faces[f].geo.cl[0]->geo.node_list[1] &&
							faces[f].geo.node_list[0]	!= faces[f].geo.cl[0]->geo.node_list[2])
						{
							faces[f].geo.cl[0]->geo.node_list[3]	= faces[f].geo.node_list[0];
							filled_nodes[faces[f].geo.cl[0]->geo.ID] = 4;

							faces[f].geo.cl[0]->geo.face_list[1]	= &faces[f];
							face_counter[faces[f].geo.cl[0]->geo.ID] =2;
						}
						else if(faces[f].geo.node_list[1]	!= faces[f].geo.cl[0]->geo.node_list[0] &&
								faces[f].geo.node_list[1]	!= faces[f].geo.cl[0]->geo.node_list[1] &&
								faces[f].geo.node_list[1]	!= faces[f].geo.cl[0]->geo.node_list[2])
						{
							faces[f].geo.cl[0]->geo.node_list[3]	= faces[f].geo.node_list[1];
							filled_nodes[faces[f].geo.cl[0]->geo.ID] = 4;

							faces[f].geo.cl[0]->geo.face_list[1]	= &faces[f];
							face_counter[faces[f].geo.cl[0]->geo.ID] =2;
						}
						else if(faces[f].geo.node_list[2]	!= faces[f].geo.cl[0]->geo.node_list[0] &&
								faces[f].geo.node_list[2]	!= faces[f].geo.cl[0]->geo.node_list[1] &&
								faces[f].geo.node_list[2]	!= faces[f].geo.cl[0]->geo.node_list[2])
						{
							faces[f].geo.cl[0]->geo.node_list[3]	= faces[f].geo.node_list[2];
							filled_nodes[faces[f].geo.cl[0]->geo.ID] = 4;

							faces[f].geo.cl[0]->geo.face_list[1]	= &faces[f];
							face_counter[faces[f].geo.cl[0]->geo.ID] = 2;
						}
					}
					else if(filled_nodes[faces[f].geo.cl[0]->geo.ID] == 4)	// SECOND FILLING
					{
						if (face_counter[faces[f].geo.cr[0]->geo.ID] != 4)
						{
							bool flag = false;
							for (int j = 0; j <= face_counter[faces[f].geo.cl[0]->geo.ID]-1; j++)
							{
								if (faces[f].geo.cl[0]->geo.face_list[j] == &faces[j])	flag = true;
							}

							if (flag != true)
							{
								int j	= face_counter[faces[f].geo.cl[0]->geo.ID];
								faces[f].geo.cl[0]->geo.face_list[j] = &faces[j];
								face_counter[faces[f].geo.cl[0]->geo.ID]	+= 1;
							}
						}
					}
				}
				else
				{
					throw ExceptionCellConstruction (FromHere(), "PROBLEM DURING CELL MESH CONSTRUCTION: Tetrajedrpm cell --ONLY TRIANGLE TYPE FACE IS ALLOWEED. PLEASE CHECK THE MESH DATA");
				}
				break;

			case CellType::quadrilateral:
				// QUADRILATERAL
				if (faces[f].geo.type == FaceType::f_line)
				{
					if (filled_nodes[faces[f].geo.cl[0]->geo.ID] == 0)		// FIRST FILLING
					{
						faces[f].geo.cl[0]->geo.node_list[0]	= faces[f].geo.node_list[0];
						faces[f].geo.cl[0]->geo.node_list[1]	= faces[f].geo.node_list[1];
						filled_nodes[faces[f].geo.cl[0]->geo.ID] = 2;

						faces[f].geo.cl[0]->geo.face_list[0]	= &faces[f];
						face_counter[faces[f].geo.cl[0]->geo.ID] = 1;
					}
					else if (filled_nodes[faces[f].geo.cl[0]->geo.ID] == 2)	// SECOND FILLING
					{
						if (faces[f].geo.node_list[0] != faces[f].geo.cl[0]->geo.node_list[0] &&
							faces[f].geo.node_list[0] != faces[f].geo.cl[0]->geo.node_list[1] &&
							faces[f].geo.node_list[1] != faces[f].geo.cl[0]->geo.node_list[0] &&
							faces[f].geo.node_list[1] != faces[f].geo.cl[0]->geo.node_list[1])
						{
							faces[f].geo.cl[0]->geo.node_list[2]	= faces[f].geo.node_list[0];
							faces[f].geo.cl[0]->geo.node_list[3]	= faces[f].geo.node_list[1];
							filled_nodes[faces[f].geo.cl[0]->geo.ID] = 4;

							faces[f].geo.cl[0]->geo.face_list[2]	= &faces[f];
							face_counter[faces[f].geo.cl[0]->geo.ID] += 1;
						}
						else if(faces[f].geo.node_list[0] == faces[f].geo.cl[0]->geo.node_list[1] &&
								faces[f].geo.node_list[1] != faces[f].geo.cl[0]->geo.node_list[2])
						{
							faces[f].geo.cl[0]->geo.face_list[1]	= &faces[f];
							face_counter[faces[f].geo.cl[0]->geo.ID] += 1;
						}
						else if(faces[f].geo.node_list[1] == faces[f].geo.cl[0]->geo.node_list[0] &&
								faces[f].geo.node_list[0] != faces[f].geo.cl[0]->geo.node_list[3])
						{
							faces[f].geo.cl[0]->geo.face_list[3]	= &faces[f];
							face_counter[faces[f].geo.cl[0]->geo.ID] += 1;
						}
					}
					else if (filled_nodes[faces[f].geo.cl[0]->geo.ID] == 4)
					{
						if (face_counter[faces[f].geo.cl[0]->geo.ID] >= 0 && face_counter[faces[f].geo.cl[0]->geo.ID] < 4)
						{
							if (faces[f].geo.node_list[0]	== faces[f].geo.cl[0]->geo.node_list[1] &&
								faces[f].geo.node_list[1]	== faces[f].geo.cl[0]->geo.node_list[2])
							{
								faces[f].geo.cl[0]->geo.face_list[1]	= &faces[f];
								face_counter[faces[f].geo.cl[0]->geo.ID] += 1;
							}
							else if(faces[f].geo.node_list[0]	== faces[f].geo.cl[0]->geo.node_list[3] &&
									faces[f].geo.node_list[1]	== faces[f].geo.cl[0]->geo.node_list[0])
							{
								faces[f].geo.cl[0]->geo.face_list[3]	= &faces[f];
								face_counter[faces[f].geo.cl[0]->geo.ID] += 1;
							}
						}
					}
				}
				else
				{
					throw ExceptionCellConstruction (FromHere(), "PROBLEM DURING CELL MESH CONSTRUCTION: QUADRILATERAL cell --ONLY Line TYPE FACE IS ALLOWEED. PLEASE CHECK THE MESH DATA");
				}
				break;

			case CellType::hexahedron:
				if (faces[f].geo.type == FaceType::f_quadrilateral)
				{
					if (filled_nodes[faces[f].geo.cl[0]->geo.ID] == 0)
					{
						faces[f].geo.cl[0]->geo.node_list[0]	= faces[f].geo.node_list[0];
						faces[f].geo.cl[0]->geo.node_list[1]	= faces[f].geo.node_list[3];
						faces[f].geo.cl[0]->geo.node_list[2]	= faces[f].geo.node_list[2];
						faces[f].geo.cl[0]->geo.node_list[3]	= faces[f].geo.node_list[1];
						filled_nodes[faces[f].geo.cl[0]->geo.ID] = 4;

						faces[f].geo.cl[0]->geo.face_list[0]	= &faces[f];
						face_counter[faces[f].geo.cl[0]->geo.ID] = 1;
					}
					else if(filled_nodes[faces[f].geo.cl[0]->geo.ID] < 8)
					{
						if (faces[f].geo.cl[0]->geo.node_list[0]	== faces[f].geo.node_list[0])
						{
							if (faces[f].geo.cl[0]->geo.node_list[1]	== faces[f].geo.node_list[1])
							{
								faces[f].geo.cl[0]->geo.node_list[4]	= faces[f].geo.node_list[2];
								faces[f].geo.cl[0]->geo.node_list[5]	= faces[f].geo.node_list[3];
								filled_nodes[faces[f].geo.cl[0]->geo.ID]			+= 2;
							}
						}
						else if (faces[f].geo.cl[0]->geo.node_list[0]	== faces[f].geo.node_list[3])
						{
							if (faces[f].geo.cl[0]->geo.node_list[1]	== faces[f].geo.node_list[0])
							{
								faces[f].geo.cl[0]->geo.node_list[4]	= faces[f].geo.node_list[2];
								faces[f].geo.cl[0]->geo.node_list[5]	= faces[f].geo.node_list[1];
								filled_nodes[faces[f].geo.cl[0]->geo.ID]			+= 2;
							}
						}
						else if (faces[f].geo.cl[0]->geo.node_list[0]	== faces[f].geo.node_list[2])
						{
							if (faces[f].geo.cl[0]->geo.node_list[1]	== faces[f].geo.node_list[3])
							{
								faces[f].geo.cl[0]->geo.node_list[4]	= faces[f].geo.node_list[1];
								faces[f].geo.cl[0]->geo.node_list[5]	= faces[f].geo.node_list[0];
								filled_nodes[faces[f].geo.cl[0]->geo.ID] 			+= 2;
							}
						}
						else if (faces[f].geo.cl[0]->geo.node_list[0]	== faces[f].geo.node_list[1])
						{
							if (faces[f].geo.cl[0]->geo.node_list[1]	== faces[f].geo.node_list[2])
							{
								faces[f].geo.cl[0]->geo.node_list[4]	= faces[f].geo.node_list[0];
								faces[f].geo.cl[0]->geo.node_list[5]	= faces[f].geo.node_list[3];
								filled_nodes[faces[f].geo.cl[0]->geo.ID]			+= 2;
							}
						}
						else if (faces[f].geo.cl[0]->geo.node_list[3] == faces[f].geo.node_list[0])
						{
							if (faces[f].geo.cl[0]->geo.node_list[2]	== faces[f].geo.node_list[3])
							{
								faces[f].geo.cl[0]->geo.node_list[6]	= faces[f].geo.node_list[2];
								faces[f].geo.cl[0]->geo.node_list[7]	= faces[f].geo.node_list[2];
								filled_nodes[faces[f].geo.cl[0]->geo.ID] 			+= 2;
							}
						}
						else if (faces[f].geo.cl[0]->geo.node_list[3]	== faces[f].geo.node_list[3])
						{
							if (faces[f].geo.cl[0]->geo.node_list[2]	== faces[f].geo.node_list[2])
							{
								faces[f].geo.cl[0]->geo.node_list[6]	= faces[f].geo.node_list[1];
								faces[f].geo.cl[0]->geo.node_list[7]	= faces[f].geo.node_list[0];
								filled_nodes[faces[f].geo.cl[0]->geo.ID]			+= 2;
							}
						}
						else if (faces[f].geo.cl[0]->geo.node_list[3]	== faces[f].geo.node_list[2])
						{
							if (faces[f].geo.cl[0]->geo.node_list[2]	== faces[f].geo.node_list[1])
							{
								faces[f].geo.cl[0]->geo.node_list[6]	= faces[f].geo.node_list[3];
								faces[f].geo.cl[0]->geo.node_list[7]	= faces[f].geo.node_list[0];
								filled_nodes[faces[f].geo.cl[0]->geo.ID]			+= 2;
							}
						}
						else if (faces[f].geo.cl[0]->geo.node_list[3] == faces[f].geo.node_list[1])
						{
							if (faces[f].geo.cl[0]->geo.node_list[2]	== faces[f].geo.node_list[0])
							{
								faces[f].geo.cl[0]->geo.node_list[6]	= faces[f].geo.node_list[3];
								faces[f].geo.cl[0]->geo.node_list[7]	= faces[f].geo.node_list[2];
								filled_nodes[faces[f].geo.cl[0]->geo.ID] 			+= 2;
							}
						}
					}
					else if(filled_nodes[faces[f].geo.cl[0]->geo.ID] > 8)
					{
						throw ExceptionCellConstruction (FromHere(), "PROBLEM DURING CELL MESH CONSTRUCTION! FOR HEXAHEDRON CELL, NUMBER OF NODE SHOULD BE 8. PLEASE CHECK THE MESH DATA!!");
					}
				}
				else
				{
					throw ExceptionCellConstruction (FromHere(), "PROBLEM DURING CELL MESH CONSTRUCTION! FOR HEXAHEDRON CELL, ONLY QUADRILATERAL TYPE FACE IS ALLOWEED. PLEASE CHECK THE MESH DATA!!");
				}
				break;

			case CellType::pyramid:
				break;

			case CellType::wedge:
				break;

			}
		}
	}



	for (int c = 1; c <= NCM; c++)	cells[c].geo.BC	== BCType::interior;

	// Calculate cell center
	for (int c = 1; c <= NCM; c++)
	{
		for (int k = 0; k <= ND-1; k++)
		{
			cells[c].geo.x[k]	= 0.0;
			for (int n = 0; n <= cells[c].geo.NN-1; n++)
			{
				cells[c].geo.x[k]	+= cells[c].geo.node_list[n]->geo.x[k];
			}

			cells[c].geo.x[k]	/= cells[c].geo.NN;


			if (cells[c].geo.x[k]	!= cells[c].geo.x[k]) throw Common::ExceptionNaNValue (FromHere(), "Nan value for cell center location ");
			if (fabs(cells[c].geo.x[k]) == numeric_limits<double>::infinity()) throw Common::ExceptionInfiniteValue (FromHere(), "Infinite Value for cell center location");
		}
	}


	// Calculate cell area/volume
	for (int c = 1; c <= NCM; c++)
	{
		switch (cells[c].geo.type)
		{
		case CellType::triangle:
			cells[c].geo.S	= Math::CalAreaTriangle(cells[c].geo.node_list[0]->geo.x, cells[c].geo.node_list[1]->geo.x, cells[c].geo.node_list[2]->geo.x);
			break;

		case CellType::tetrahedron:
			cells[c].geo.S	= Math::CalVolumeTetrahedron(cells[c].geo.node_list[0]->geo.x, cells[c].geo.node_list[1]->geo.x, cells[c].geo.node_list[2]->geo.x, cells[c].geo.node_list[3]->geo.x);
			break;

		case CellType::quadrilateral:
			cells[c].geo.S	= Math::CalAreaQuadrilateral(cells[c].geo.node_list[0]->geo.x, cells[c].geo.node_list[1]->geo.x, cells[c].geo.node_list[2]->geo.x, cells[c].geo.node_list[3]->geo.x);
			break;

		case CellType::hexahedron:
			cells[c].geo.S	= Math::CalVolumeHexahedron(cells[c].geo.node_list[0]->geo.x, cells[c].geo.node_list[1]->geo.x, cells[c].geo.node_list[2]->geo.x, cells[c].geo.node_list[3]->geo.x, cells[c].geo.node_list[4]->geo.x, cells[c].geo.node_list[5]->geo.x, cells[c].geo.node_list[6]->geo.x, cells[c].geo.node_list[7]->geo.x);
			break;

		case CellType::pyramid:
			cells[c].geo.S	= Math::CalVolumePyramid(cells[c].geo.node_list[0]->geo.x, cells[c].geo.node_list[1]->geo.x, cells[c].geo.node_list[2]->geo.x, cells[c].geo.node_list[3]->geo.x, cells[c].geo.node_list[4]->geo.x);
			break;

		case CellType::wedge:
			cells[c].geo.S	= Math::CalVolumeWedge(cells[c].geo.node_list[0]->geo.x, cells[c].geo.node_list[1]->geo.x, cells[c].geo.node_list[2]->geo.x, cells[c].geo.node_list[3]->geo.x, cells[c].geo.node_list[4]->geo.x, cells[c].geo.node_list[4]->geo.x);
			break;

		}
	}


	// Calculate characteristic length
	for (int c = 1; c <= NCM; c++)
	{
		vector <double> Xf1(ND, 0.0);
		vector <double> Xf2(ND, 0.0);
		vector <double> Xf3(ND, 0.0);
		vector <double> Xf4(ND, 0.0);

		double aux;
		double length;

		switch (cells[c].geo.type)
		{
		case CellType::triangle:
			cells[c].geo.characteristic_length	= sqrt(cells[c].geo.S);
			break;

		case CellType::tetrahedron:
			for (int k = 0; k <= ND-1; k++)	Xf1[k]	= (cells[c].geo.node_list[0]->geo.x[k] + cells[c].geo.node_list[1]->geo.x[k] + cells[c].geo.node_list[2]->geo.x[k]) / 3.0;
			for (int k = 0; k <= ND-1; k++)	Xf2[k]	= (cells[c].geo.node_list[0]->geo.x[k] + cells[c].geo.node_list[1]->geo.x[k] + cells[c].geo.node_list[3]->geo.x[k]) / 3.0;
			for (int k = 0; k <= ND-1; k++)	Xf3[k]	= (cells[c].geo.node_list[1]->geo.x[k] + cells[c].geo.node_list[2]->geo.x[k] + cells[c].geo.node_list[3]->geo.x[k]) / 3.0;
			for (int k = 0; k <= ND-1; k++)	Xf4[k]	= (cells[c].geo.node_list[0]->geo.x[k] + cells[c].geo.node_list[3]->geo.x[k] + cells[c].geo.node_list[2]->geo.x[k]) / 3.0;

			aux	= 0.0;
			for (int k = 0; k <= ND-1; k++)	aux += pow(cells[c].geo.x[k] - Xf1[k], 2.0);
			length	= sqrt(aux);

			aux = 0.0;
			for (int k = 0; k <= ND-1; k++)	aux += pow(cells[c].geo.x[k] - Xf2[k], 2.0);
			aux 	= sqrt(aux);
			length	= Math::fmin<double>(length, aux);

			aux = 0.0;
			for (int k = 0; k <= ND-1; k++)	aux += pow(cells[c].geo.x[k] - Xf3[k], 2.0);
			aux 	= sqrt(aux);
			length	= Math::fmin<double>(length, aux);


			aux = 0.0;
			for (int k = 0; k <= ND-1; k++)	aux += pow(cells[c].geo.x[k] - Xf4[k], 2.0);
			aux 	= sqrt(aux);
			length	= Math::fmin<double>(length, aux);

			cells[c].geo.characteristic_length = 2.0* length;
			break;

		case CellType::quadrilateral:
			for (int k = 0; k <= ND-1; k++)	Xf1[k]	= cells[c].geo.node_list[1]->geo.x[k] - cells[c].geo.node_list[0]->geo.x[k];
			for (int k = 0; k <= ND-1; k++)	Xf2[k]	= cells[c].geo.node_list[3]->geo.x[k] - cells[c].geo.node_list[0]->geo.x[k];

			length	= 0.0;	aux		= 0.0;
			for (int k = 0; k <= ND-1; k++)
			{
				length	+= pow(Xf1[k], 2.0);
				aux		+= pow(Xf2[k], 2.0);
			}

			length	= Math::fmin<double>(length, aux);
			cells[c].geo.characteristic_length = sqrt(length);
			break;


		case CellType::hexahedron:
			for (int k = 0; k <= ND-1; k++)	Xf1[k]	= cells[c].geo.node_list[1]->geo.x[k] - cells[c].geo.node_list[0]->geo.x[k];
			for (int k = 0; k <= ND-1; k++)	Xf2[k]	= cells[c].geo.node_list[2]->geo.x[k] - cells[c].geo.node_list[0]->geo.x[k];
			for (int k = 0; k <= ND-1; k++)	Xf3[k]	= cells[c].geo.node_list[3]->geo.x[k] - cells[c].geo.node_list[0]->geo.x[k];

			length	= 0.0;	aux		= 0.0;
			for (int k = 0; k <= ND-1; k++)
			{
				length	+= pow(Xf1[k], 2.0);
				aux		+= pow(Xf2[k], 2.0);
			}
			length	= Math::fmin<double>(length, aux);

			aux = 0.0;
			for (int k = 0; k <= ND-1; k++)	aux += pow(Xf3[k], 2.0);
			length	= Math::fmin<double>(length, aux);

			cells[c].geo.characteristic_length = sqrt(length);
			break;

		case CellType::pyramid:
			length = pow(cells[c].geo.S, 1.0/3.0);
			cells[c].geo.characteristic_length = sqrt(length);
			break;

		case CellType::wedge:
			for (int k = 0; k <= ND-1; k++)	Xf1[k]	= cells[c].geo.node_list[1]->geo.x[k] - cells[c].geo.node_list[0]->geo.x[k];
			for (int k = 0; k <= ND-1; k++)	Xf2[k]	= cells[c].geo.node_list[2]->geo.x[k] - cells[c].geo.node_list[0]->geo.x[k];
			for (int k = 0; k <= ND-1; k++)	Xf3[k]	= cells[c].geo.node_list[3]->geo.x[k] - cells[c].geo.node_list[0]->geo.x[k];

			aux = 	pow(0.5*(Xf1[1]*Xf2[2] - Xf1[2]*Xf2[1]), 2.0);
			aux	+= 	pow(-0.5*(Xf1[0]*Xf2[2] - Xf1[2]*Xf2[0]), 2.0);
			aux	+=	pow(0.5*(Xf1[0]*Xf2[1] - Xf1[1]*Xf2[0]), 2.0);
			length	= pow(aux, 0.25);

			aux = 0.0;
			for (int k = 0; k <= ND-1; k++)	aux += pow(Xf3[k], 2.0);
			aux = sqrt(aux);

			length	= Math::fmin<double>(length, aux);
			cells[c].geo.characteristic_length = length;
			break;
		}
	}
}





// Create/processing Ghost-cell data
void Grid::processingGhostData()
{
	int counter = 0;

	for (int f = 1; f <= NFM; f++)
	{
		if (faces[f].geo.cr[0]	== NULL)
		{
			counter++;
			faces[f].geo.cr[0]	= &cells_ghost[counter];

			cells_ghost[counter].geo.ID	= -counter;
			cells_ghost[counter].geo.allocate(ND, CellType::ghost);

			cells_ghost[counter].geo.NN	= faces[f].geo.NN;
			cells_ghost[counter].geo.node_list.resize(faces[f].geo.NN);
			for (int n = 0; n <= faces[f].geo.NN-1; n++)	cells_ghost[counter].geo.node_list[n]	= faces[f].geo.node_list[n];

			for (int k = 0; k <= ND-1; k++)	cells_ghost[counter].geo.x[k]	= 2.0*faces[f].geo.x[k] - faces[f].geo.cl[0]->geo.x[k];
			cells_ghost[counter].geo.S						= faces[f].geo.S;
			cells_ghost[counter].geo.characteristic_length	= faces[f].geo.S;

			cells_ghost[counter].geo.face_list[0]	= &faces[f];
			cells_ghost[counter].geo.BC				= faces[f].geo.BC;
			cells_ghost[counter].geo.type			= CellType::ghost;
		}
		else if (faces[f].geo.cl[0]	== NULL)
		{
			counter++;
			faces[f].geo.cl[0]	= faces[f].geo.cr[0];
			faces[f].geo.cr[0]	= &cells_ghost[counter];

			for (int k1 = 0; k1 <= ND-1; k1++)
				for (int k2 = 0; k2 <= ND-1; k2++)
					faces[f].geo.n[k1][k2]	= -faces[f].geo.n[k1][k2];

			cells_ghost[counter].geo.ID	= -counter;
			cells_ghost[counter].geo.allocate(ND, CellType::ghost);

			cells_ghost[counter].geo.NN	= faces[f].geo.NN;
			cells_ghost[counter].geo.node_list.resize(faces[f].geo.NN);
			for (int n = 0; n <= faces[f].geo.NN-1; n++)	cells_ghost[counter].geo.node_list[n]	= faces[f].geo.node_list[n];

			for (int k = 0; k <= ND-1; k++)	cells_ghost[counter].geo.x[k]	= 2.0*faces[f].geo.x[k] - faces[f].geo.cl[0]->geo.x[k];
			cells_ghost[counter].geo.S						= faces[f].geo.S;
			cells_ghost[counter].geo.characteristic_length	= faces[f].geo.S;

			cells_ghost[counter].geo.face_list[0]	= &faces[f];
			cells_ghost[counter].geo.BC				= faces[f].geo.BC;
			cells_ghost[counter].geo.type			= CellType::ghost;

		}
	}


	if (counter != NGM)
	{
		throw ExceptionGridDataMismatch (FromHere(), "PROBLEM IN the creation of ghost-Cell DATA. TOTAL NUMBER OF created ghost-Cell DATA DOES NOT MATHCH WITH MESH INFOMATION DATA");
	}
}






// Grid processing function
void Grid::processingGridData(const double mesh_factor, bool is_axisymmetric, bool extended_stencill)
{
	processingNodeData(mesh_factor, is_axisymmetric);
	processingFaceData();
	processingCellData();
	processingGhostData();

	find_stencil(extended_stencill);
	calculate_distance_to_wall();
}



}
}
