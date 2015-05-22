/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_utilities_node_value.cpp
 * 			-  
 *  
 */






#include <mkl.h>
#include <vector>
#include <limits>
#include "omp.h"

#include "../../GRID/include/OP2A_grid.hpp"
#include "../../MATRIX/include/CMatrix.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../../DATA/include/OP2A_data.hpp"


using namespace std;

void calculate_node_value_find_make_lists(GRID_CLASS &grid, vector < vector<int> >	&node_shared_cell_list, vector < vector<double> >	&node_shared_cell_weighting, int nt)
{
	kmp_set_defaults("KMP_AFFINITY = scatter");

	// Initialize values
#pragma omp parallel for num_threads(nt)
	for (int n = 0; n <= grid.NNM; n++)	node_shared_cell_list[n][0]	= 0;


	// Count number of shared cells per each node
	for (int f = 0; f <= grid.NFM-1; f++)
	{
		int cl	= grid.faces.data_ptr[f]->cl[0];
		int cr	= grid.faces.data_ptr[f]->cr[0];

		for (int k = 0; k <= grid.faces.data_ptr[f]->NN-1; k++)
		{
			int flag;
			int n		= grid.faces.data_ptr[f]->node[k];
			int index	= node_shared_cell_list[n][0];

			flag = 0;
			for (int i = 1; i <= node_shared_cell_list[n][0]; i++)	if (node_shared_cell_list[n][i]	== cl)	flag = 1;
			if (flag == 0)
			{
				index++;
				node_shared_cell_list[n][index]	= cl;
			}


			flag = 0;
			for (int i = 1; i <= node_shared_cell_list[n][0]; i++)	if (node_shared_cell_list[n][i]	== cr)	flag = 1;
			if (flag == 0)
			{
				index++;
				node_shared_cell_list[n][index]	= cr;
			}

			node_shared_cell_list[n][0]	= index;
		}
	}


	// Use unit weighting value
#pragma omp parallel for num_threads(nt)
	for (int n = 1; n <= grid.NNM; n++)
	{
		for (int c = 1; c <= node_shared_cell_list[n][0]; c++)
		{
			node_shared_cell_weighting[n][c] = 1.0 / node_shared_cell_list[n][0];
		}
	}
}



void calculate_node_value(GRID_CLASS &grid, SOL_CLASS_DATA &cell_Q,  SOL_CLASS_DATA &cellg_Q,  SOL_CLASS_DATA &node_Q,
						vector < vector<int> >	&node_shared_cell_list, vector < vector<double> >	&node_shared_cell_weighting,
						vector<int> variable_index, int num_variable, int nt)
{
	kmp_set_defaults("KMP_AFFINITY = scatter");

#pragma omp parallel for num_threads(nt)
	for (int n = 1; n <= grid.NNM; n++)
	{
		int n_ptr	= grid.nodes.whereis[n];

		for (int v = 0; v <= num_variable-1; v++)
		{
			int j	= variable_index[v];
			node_Q.data_ptr[n_ptr]->data[v]	= 0.0;

			for (int i = 1; i <= node_shared_cell_list[n][0]; i++)
			{
				int c_ptr;
				int c = node_shared_cell_list[n][i];

				if (c > 0)
				{
					c_ptr	= grid.cells.whereis[c];
					node_Q.data_ptr[n_ptr]->data[v]	+= cell_Q.data_ptr[c_ptr]->data[j]*node_shared_cell_weighting[n][i];
				}
				else
				{
					c_ptr	= grid.cells_ghost.whereis[-c];
					node_Q.data_ptr[n_ptr]->data[v]	+= cellg_Q.data_ptr[c_ptr]->data[j]*node_shared_cell_weighting[n][i];
				}
			}
		}
	}
}
