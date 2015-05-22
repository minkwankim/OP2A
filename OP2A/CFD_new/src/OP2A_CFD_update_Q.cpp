/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 16, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_update_Q.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_CFD_Time_integration.hpp"

void CFD_update_Q(GRID_CLASS &grid, SOL_CFD &Solution, int nt)
{

	int ne	= Solution.setup.NS + Solution.setup.ND;

#pragma omp parallel for num_threads(nt)
	for (int c = 0; c <= grid.NCM-1; c++)
	{
		for (int i = 0;  i<= Solution.setup.VAR-1; i++)
			Solution.Q_new[c][i]	= Solution.Qc[c][i] + Solution.dQ[c][i];

		// Preliminary Error Check
		//	A. Density
		for (int s = 0; s <= Solution.setup.NS-1; s++)
		{
			if (Solution.Q_new[c][s] < 0.0 && fabs(Solution.Q_new[c][s]) <= 1.0e-15)
			{
				Solution.Q_new[c][s] = 0.0;
			}

			if (error_check_double_pos(Solution.Q_new[c][s]) == true)
			{
				Error_message_type	error_message;
				error_message.module_name	="CFD VARIABLE UPDATE";
				error_message.location_primary_name	= "Cell ID";
				error_message.location_primary		= grid.cells.data_ptr[c]->ID;
				error_message.location_secondary_name	= "VARIABLE";
				error_message.location_secondary		= s;
				error_message.message	= " Negative Density!";
				error_message.print_message();
			}

		}

		// B. Energy
		for (int m = 0; m <= Solution.setup.NE-1; m++)
		{
			if (error_check_double_pos(Solution.Q_new[c][ne+m]) == true)
			{
				Error_message_type	error_message;
				error_message.module_name	="CFD VARIABLE UPDATE";
				error_message.location_primary_name	= "Cell ID";
				error_message.location_primary		= grid.cells.data_ptr[c]->ID;
				error_message.location_secondary_name	= "VARIABLE";
				error_message.location_secondary		= m;
				error_message.message	= " Negative Energy!";
				error_message.print_message();
			}
		}
	}
}

