/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 9, 2015
 *      			Author: Minkwan Kim
 *
 * GridStencilFinder.cpp
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
#include "Common/include/MultiDimension.hpp"


#include "Math/include/AreaCalculation.hpp"
#include "Math/include/MathMisc.hpp"
#include "Math/include/OP2A_Vector.hpp"



namespace OP2A{
namespace GRID{

void Grid::find_stencil(bool extended_stencill)
{

	// Step 1: Assign Basic variables
	for (int n = 1; n <= NNM; n++)	nodes[n].geo.NSC = 0;


	// Step 2: FIND OUT EACH NODES ARE BELONG TO HOW MANY CELLS
	for (int c = 1; c <= NCM; c++)
	{
		for (int n = 0; n <= cells[c].geo.NN-1; n++)	cells[c].geo.node_list[n]->geo.NSC++;
	}

	for (int c = 1; c <= NGM; c++)
	{
		for (int n = 0; n <= cells_ghost[c].geo.NN-1; n++)	cells_ghost[c].geo.node_list[n]->geo.NSC++;
	}

	int max_cell_sharing_node = 0;
	for (int n = 1; n <= NNM; n++)
	{
		nodes[n].geo.WeightShairedCell.resize(nodes[n].geo.NSC);
		nodes[n].geo.CellList.resize(nodes[n].geo.NSC);

		max_cell_sharing_node = Math::fmax<int>(nodes[n].geo.NSC, max_cell_sharing_node);
	}

	for (int n = 1; n <= NNM; n++)	nodes[n].geo.NSC = 0;



	// Step 3: GET SHARED CELL ID NUMBER FOR EACH NODES
	for (int c = 1; c <= NCM; c++)
	{
		for (int n = 0; n <= cells[c].geo.NN-1; n++)
		{
			int index = cells[c].geo.node_list[n]->geo.NSC;

			cells[c].geo.node_list[n]->geo.WeightShairedCell[index]	= 1.0;
			cells[c].geo.node_list[n]->geo.CellList[index]	= &cells[c];
			cells[c].geo.node_list[n]->geo.NSC++;
		}
	}

	for (int c = 1; c <= NGM; c++)
	{
		for (int n = 0; n <= cells_ghost[c].geo.NN-1; n++)
		{
			int index = cells_ghost[c].geo.node_list[n]->geo.NSC;

			cells_ghost[c].geo.node_list[n]->geo.WeightShairedCell[index]	= 1.0;
			cells_ghost[c].geo.node_list[n]->geo.CellList[index]			= &cells_ghost[c];
			cells_ghost[c].geo.node_list[n]->geo.NSC++;
		}
	}

	for (int n = 1; n <= NNM; n++)
	{

		double S_temp	= 0.0;
		for (int ii = 0; ii <= nodes[n].geo.NSC-1; ii++)
		{
			double area_temp;
			if (nodes[n].geo.CellList[ii]->geo.type != CellType::ghost)
			{
				area_temp = nodes[n].geo.CellList[ii]->geo.S;
			}
			else
			{
				area_temp = nodes[n].geo.CellList[ii]->geo.face_list[0]->geo.cl[0]->geo.S;
			}

			nodes[n].geo.WeightShairedCell[ii]	= area_temp;
			S_temp += area_temp;
		}

		for (int ii = 0; ii <= nodes[n].geo.NSC-1; ii++)	nodes[n].geo.WeightShairedCell[ii] /= S_temp;
	}



	// Step 4: GET NEIGHBOR CELL List
	for (int c = 1; c <= NCM; c++)
	{
		if (cells[c].geo.neighbor_list.size() != cells[c].geo.NF)	cells[c].geo.neighbor_list.resize(cells[c].geo.NF);

		for (int f = 0; f <= cells[c].geo.NF-1; f++)
		{
			if (cells[c].geo.face_list[f]->geo.cl[0] == &cells[c])	cells[c].geo.neighbor_list[f]	= cells[c].geo.face_list[f]->geo.cr[0];
			else													cells[c].geo.neighbor_list[f]	= cells[c].geo.face_list[f]->geo.cl[0];
		}
	}





	// Step 5: Find Stencil
	double dot;
	double dot2;

	Cell *trial_cell;
	Cell *left_cell;
	Cell *right_cell;

	for (int f = 1; f <= NFM; f++)
	{
		left_cell	= faces[f].geo.cl[0];
		right_cell	= faces[f].geo.cr[0];


		// Step 5.1: FIND CLL / CL1 / CL2
		double max_dot = -1.0;
		double max_dot2	= -1.0;
		double min_dot2	= 1.0;


		if (extended_stencill == false)
		{
			// GET NORMAL/TANGENTIAL VECTOR TOWARD TO CL
			Math::VECTOR normal(faces[f].geo.n[0]);		normal *= -1.0;

			// CALCULATE DISTANCE TO CELL-CENTER OF EACH NEIGHBORING CELL
			for (int ff = 0; ff <= left_cell->geo.NF-1; ff++)
			{
				trial_cell	= left_cell->geo.neighbor_list[ff];

				Math::VECTOR	Xfc(left_cell->geo.x, trial_cell->geo.x);	Xfc.normalize();

				dot 	= Math::VectorDotProduct(normal, Xfc);

				if (dot > max_dot)
				{
					max_dot = dot;
					faces[f].geo.cl[StencilLabel::CLL]	= trial_cell;
				}
			}

			if (faces[f].geo.cl[StencilLabel::CLL] == NULL)
			{
				faces[f].geo.cl[StencilLabel::CLL] = faces[f].geo.cl[StencilLabel::CL];
			}
		}
		else
		{
			// GET NORMAL/TANGENTIAL VECTOR TOWARD TO CL
			Math::VECTOR normal(faces[f].geo.n[0]);		normal *= -1.0;
			Math::VECTOR tangential(faces[f].geo.n[0]);

			for (int ff = 0; ff <= left_cell->geo.NF-1; ff++)
			{
				trial_cell	= left_cell->geo.neighbor_list[ff];

				Math::VECTOR	Xfc(left_cell->geo.x, trial_cell->geo.x);	Xfc.normalize();


				dot 	= Math::VectorDotProduct(normal, Xfc);
				dot2	= Math::VectorDotProduct(tangential, Xfc);

				if (dot > max_dot)
				{
					max_dot = dot;
					faces[f].geo.cl[StencilLabel::CLL]	= trial_cell;
				}

				if (dot2 > max_dot2)
				{
					max_dot2	= dot2;
					faces[f].geo.cl[StencilLabel::CLupper]	= trial_cell;
				}

				if (dot2 < min_dot2)
				{
					min_dot2	= dot2;
					faces[f].geo.cl[StencilLabel::CLlower]	= trial_cell;
				}
			}


			if (faces[f].geo.cl[StencilLabel::CLL] == NULL)
			{
				faces[f].geo.cl[StencilLabel::CLL] = faces[f].geo.cl[StencilLabel::CL];
			}

			if (faces[f].geo.cl[StencilLabel::CLupper] == NULL)
			{
				faces[f].geo.cl[StencilLabel::CLupper] = faces[f].geo.cl[StencilLabel::CL];
			}

			if (faces[f].geo.cl[StencilLabel::CLlower] == NULL)
			{
				faces[f].geo.cl[StencilLabel::CLlower] = faces[f].geo.cl[StencilLabel::CL];
			}
		}


		// Step 5.2: FIND CRR
		if (right_cell->geo.type != CellType::ghost)	// NOT A GHOST CELL
		{
			max_dot 	= -1.0;
			max_dot2	= -1.0;
			min_dot2	= 1.0;

			if (extended_stencill == false)
			{
				// GET NORMAL/TANGENTIAL VECTOR TOWARD TO CL
				Math::VECTOR normal(faces[f].geo.n[0]);

				// CALCULATE DISTANCE TO CELL-CENTER OF EACH NEIGHBORING CELL
				for (int ff = 0; ff <= right_cell->geo.NF-1; ff++)
				{
					trial_cell	= right_cell->geo.neighbor_list[ff];

					Math::VECTOR	Xfc(right_cell->geo.x, trial_cell->geo.x);	Xfc.normalize();

					dot 	= Math::VectorDotProduct(normal, Xfc);

					if (dot > max_dot)
					{
						max_dot = dot;
						faces[f].geo.cr[StencilLabel::CRR]	= trial_cell;
					}
				}

				if (faces[f].geo.cl[StencilLabel::CRR] == NULL)
				{
					faces[f].geo.cl[StencilLabel::CRR] = faces[f].geo.cl[StencilLabel::CR];
				}
			}
			else
			{
				// GET NORMAL/TANGENTIAL VECTOR TOWARD TO CL
				Math::VECTOR normal(faces[f].geo.n[0]);
				Math::VECTOR tangential(faces[f].geo.n[0]);

				for (int ff = 0; ff <= right_cell->geo.NF-1; ff++)
				{
					trial_cell	= right_cell->geo.neighbor_list[ff];

					Math::VECTOR	Xfc(right_cell->geo.x, trial_cell->geo.x);	Xfc.normalize();


					dot 	= Math::VectorDotProduct(normal, Xfc);
					dot2	= Math::VectorDotProduct(tangential, Xfc);

					if (dot > max_dot)
					{
						max_dot = dot;
						faces[f].geo.cr[StencilLabel::CRR]	= trial_cell;
					}

					if (dot2 > max_dot2)
					{
						max_dot2	= dot2;
						faces[f].geo.cr[StencilLabel::CRupper]	= trial_cell;
					}

					if (dot2 < min_dot2)
					{
						min_dot2	= dot2;
						faces[f].geo.cr[StencilLabel::CRlower]	= trial_cell;
					}
				}

			}
		}
		else
		{
			if (extended_stencill == false)
			{
				faces[f].geo.cr[StencilLabel::CRR]	= faces[f].geo.cr[StencilLabel::CR];
			}
			else
			{
				faces[f].geo.cr[StencilLabel::CRupper]	= faces[f].geo.cr[StencilLabel::CR];
				faces[f].geo.cr[StencilLabel::CRlower]	= faces[f].geo.cr[StencilLabel::CR];
			}
		}
	}
}


}
}
