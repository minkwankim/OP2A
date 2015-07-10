/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 9, 2015
 *      			Author: Minkwan Kim
 *
 * GridDistanceWall.cpp
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
#include "../include/Exception_GridStencil.hpp"


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

void Grid::calculate_distance_to_wall()
{
	// 1. Initizlize variables
	vector<int>		walk_bound(NCM+1, 0);
	vector<int>		near_wall(NCM+1, 0);	// Near wall face ID
	for (int c = 1; c <= NCM;	c++) cells[c].geo.dist_wall	= 1.E10;

	int	n_wall_faces	= 0;
	int	n_walk_cells	= 0;


	// 2. Calculate distance for 1st layer
	for (int c = 1; c <= NGM; c++)
	{
		if (is_wall_typeBC(cells_ghost[c].geo.BC))
		{
			n_wall_faces++;

			if (cells_ghost[c].geo.face_list[0]->geo.cl[0]->geo.type == CellType::ghost)	throw ExceptionGridStencil (FromHere(), "Problem in the stencil of face!");

			walk_bound[cells_ghost[c].geo.face_list[0]->geo.cl[0]->geo.ID]	= 1;
			near_wall[cells_ghost[c].geo.face_list[0]->geo.cl[0]->geo.ID]		= cells_ghost[c].geo.face_list[0]->geo.ID;
			n_walk_cells++;

			Math::VECTOR distance_to_wall_initial(cells_ghost[c].geo.face_list[0]->geo.x, cells_ghost[c].geo.face_list[0]->geo.cl[0]->geo.x);

			cells_ghost[c].geo.face_list[0]->geo.cl[0]->geo.dist_wall = distance_to_wall_initial.length();
		}
	}

	// 2. Calculate path
	int n_pass = 0;
	int wall;
	while (n_walk_cells > 0 && n_pass <= GRID_MAX_NUM_PATH)
	{
		n_pass++;
        n_walk_cells	= 0;

        for (int c = 1; c <= NCM; c++)
    	{

        	if (walk_bound[c] == n_pass)
        	{
        		wall	= near_wall[c];	// Near wall face ID

        		for (int f = 0; f <= cells[c].geo.NF-1; f++)
        		{
        			Cell *left_cell;
        			Cell *right_cell;

        			left_cell	= cells[c].geo.face_list[f]->geo.cl[0];
        			right_cell	= cells[c].geo.face_list[f]->geo.cr[0];

        			if (&cells[c] == left_cell)	left_cell = right_cell;

        			if (left_cell->geo.type	!= CellType::ghost)
        			{
        				if (walk_bound[left_cell->geo.ID]	!= 1)
        				{
        					if (near_wall[left_cell->geo.ID]	== 0)
        					{
        						near_wall[left_cell->geo.ID]	= wall;
        						walk_bound[left_cell->geo.ID]	= n_pass + 1;
        						n_walk_cells++;


        						Math::VECTOR dist(faces[wall].geo.x, left_cell->geo.x);
        						left_cell->geo.dist_wall = dist.length();
        					}
        					else if (near_wall[left_cell->geo.ID] != wall)
        					{
        						Math::VECTOR dist(faces[wall].geo.x, left_cell->geo.x);

  								if (dist.length() < left_cell->geo.dist_wall)
								{
									near_wall[left_cell->geo.ID]		= wall;
									walk_bound[left_cell->geo.ID]	= n_pass+1;
									n_walk_cells++;

									left_cell->geo.dist_wall	= dist.length();
								}
        					}
        				}
        			}
        		}
        	}
    	}
	}



	for (int c = 1; c <= NGM-1; c++)
	{
		cells_ghost[c].geo.dist_wall	= -cells_ghost[c].geo.face_list[0]->geo.cl[0]->geo.dist_wall;
	}


	for (int f = 1; f <= NFM; f++)
	{
		faces[f].geo.dist_wall	= 0.5 * (faces[f].geo.cl[0]->geo.dist_wall + faces[f].geo.cr[0]->geo.dist_wall);

		double d = 0.0;
        for (int k = 0; k <= ND-1; k++)	d += faces[f].geo.n[0][k]*faces[near_wall[faces[f].geo.cl[0]->geo.ID]].geo.n[0][k];
        faces[f].geo.n_dot_wall	= 1.0 - fabs(d);
	}
}





}
}


