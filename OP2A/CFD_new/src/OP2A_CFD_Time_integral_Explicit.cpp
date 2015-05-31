/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 16, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Time_integral_Explicit.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_CFD_Variable_change.hpp"
#include "../include/OP2A_CFD_reconstruct.hpp"
#include "../include/OP2A_CFD_Flux.hpp"


#include "../../UTIL/include/OP2A_utilities.hpp"
#include "Problem/include/OP2A_setup_problem.hpp"
#include "../../MATRIX/include/OPPA_Matrix.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../GRID/include/OP2A_grid.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"

void CFD_time_integration_Explicit(GRID_CLASS &grid, SOL_CFD &Solution, double dt, bool is_axis, int nt)
{
	int 	c;
	double	Vol, dt_Vol;

#pragma omp parallel for private(Vol, dt_Vol) num_threads(nt)
	for (c = 0; c <= grid.NCM-1; c++)
	{
		if (is_axis	== true)	Vol		= grid.cells.data_ptr[c]->S	* fabs(grid.cells.data_ptr[c]->x[1]);
		else					Vol		= grid.cells.data_ptr[c]->S;
		dt_Vol	= dt / Vol;


		for (int i = 0; i <= Solution.setup.VAR-1; i++)
		{
			Solution.dQ[c][i]	= -dt_Vol * Solution.Rn[c][i];

			if  (error_check_double(Solution.dQ[c][i]) == true)
			{
				Error_message_type	error_message;
				error_message.module_name	="CFD TIME INTEGRATION EXPLICIT";
				error_message.location_primary_name	= "Cell";
				error_message.location_primary		= grid.cells.data_ptr[c]->ID;
				error_message.location_secondary_name	= "VARIABLE";
				error_message.location_secondary		= i;
				error_message.message	= " Problem in time integration. Need to check residues!";
				error_message.print_message();
			}
		}
	}
}

