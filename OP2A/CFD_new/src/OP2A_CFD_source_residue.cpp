/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 30, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_source_residue.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_CFD_Variable_change.hpp"
#include "../include/OP2A_CFD_reconstruct.hpp"
#include "../include/OP2A_CFD_Flux.hpp"


#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../../SETUP/include/OP2A_setup_problem.hpp"
#include "../../MATRIX/include/OPPA_Matrix.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../GRID/include/OP2A_grid.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"

void CFD_residue_source(GRID_CLASS &grid_data, SOL_CFD &Solutions, bool is_axis, int nt)
{
#pragma omp parallel for num_threads(nt)
	for (int c = 0; c <= grid_data.NCM-1; c++)
	{
		double S = 0.0;

		if (is_axis	== true)	S	= grid_data.cells.data_ptr[c]->S * fabs(grid_data.cells.data_ptr[c]->x[1]);
		else					S	= grid_data.cells.data_ptr[c]->S;

		for (int s = 0; s <= Solutions.setup.VAR; s++)
		{
			double Source_term	= Solutions.S_source[c][s];
			Source_term	= Source_term * S;

			Solutions.Rn[c][s]	-= Source_term;
		}
	}
}
