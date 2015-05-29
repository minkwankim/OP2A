/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 15, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Flux_Inviscid_residue.cpp
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


void CFD_residue_inviscid_ver2(GRID_CLASS &grid_data, SOL_CFD &Solution, bool is_axis, int nt)
{
	// Step 1: Initialize Residue term
#pragma omp parallel for num_threads(nt)
	for (int c = 0; c <= grid_data.NCM-1; c++)
	{
		for (int i = 0; i <= Solution.setup.VAR-1; i++)	Solution.Rn[c][i]	= 0.0;
	}


	// Step 2: Calculate residues
	for (int f = 0; f <= grid_data.NFM-1; f++)
	{
		int cl = grid_data.faces.data_ptr[f]->cl[0];
		int cr = grid_data.faces.data_ptr[f]->cr[0];

		double S = 0.0;
		if (is_axis	== true)	S	= grid_data.faces.data_ptr[f]->S * fabs(grid_data.faces.data_ptr[f]->x[1]);
		else					S	= grid_data.faces.data_ptr[f]->S;

#pragma omp parallel for
		for (int k = 0; k <= Solution.setup.VAR-1; k++)
		{
			double flux	= Solution.Flux_f_inviscid[f][k] * S;

			if (cl > 0)
			{
				int cl_ptr	= grid_data.cells.whereis[cl];
				Solution.Rn[cl_ptr][k]	+= flux;
			}

			if (cr > 0)
			{
				int cr_ptr	= grid_data.cells.whereis[cr];
				Solution.Rn[cr_ptr][k]	-= flux;
			}
		}
	}


	// Step 3: Update for Axisymmetric case
	if (is_axis == true)
	{
		int ne	= Solution.setup.NS+Solution.setup.ND;

#pragma omp parallel for num_threads(nt)
		for (int c = 0;  c <= grid_data.NCM-1; c++)
		{
			double	S	= grid_data.cells.data_ptr[c]->S;
			Solution.Rn[c][Solution.setup.NS+1]	-= S * Solution.Wc[c][ne] * fsgn(grid_data.cells.data_ptr[c]->x[1]);
		}
	}


	// Check error
#pragma omp parallel for num_threads(nt)
	for (int c = 0; c <= grid_data.NCM-1; c++)
	{
		for (int i = 0; i <= Solution.setup.VAR-1; i++)
		{
			if (error_check_double(Solution.Rn[c][i]))
			{
				Error_message_type	error_message;
				error_message.module_name	="CFD-Inviscid Flux Residue";
				error_message.location_primary_name	= "Cell";
				error_message.location_primary		= c;
				error_message.location_secondary_name	= "Variable";
				error_message.location_secondary		= i;
				error_message.message	= " Problem in inviscid flux residue at cell!";
				error_message.print_message();
			}
		}
	}
}



void CFD_residue_inviscid_ver3(GRID_CLASS &grid_data, SOL_CFD &Solution, bool is_axis, int nt)
{
	// Step 1: Initialize Residue term
#pragma omp parallel for num_threads(nt)
	for (int c = 0; c <= grid_data.NCM-1; c++)
	{
		for (int i = 0; i <= Solution.setup.VAR-1; i++)	Solution.Rn[c][i]	= 0.0;
	}


	// Step 2: Calculate residues
#pragma omp parallel for num_threads(nt)
	for (int c = 0; c <= grid_data.NCM-1; c++)
	{
		for (int i = 0; i <= grid_data.cells.data_ptr[c]->NF-1; i++)
		{
			int f_ID	= grid_data.cells.data_ptr[c]->face[i];
			int f		= grid_data.faces.whereis[f_ID];

			int cl = grid_data.faces.data_ptr[f]->cl[0];
			int cr = grid_data.faces.data_ptr[f]->cr[0];

			double S = 0.0;
			if (is_axis	== true)	S	= grid_data.faces.data_ptr[f]->S * fabs(grid_data.faces.data_ptr[f]->x[1]);
			else					S	= grid_data.faces.data_ptr[f]->S;

			for (int k = 0; k <= Solution.setup.VAR-1; k++)
			{
				double flux	= Solution.Flux_f_inviscid[f][k] * S;

				if (cr == grid_data.cells.data_ptr[c]->ID)		Solution.Rn[c][k]	-= flux;
				else if (cl == grid_data.cells.data_ptr[c]->ID)	Solution.Rn[c][k]	+= flux;
			}
		}
	}


	// Step 3: Update for Axisymmetric case
	if (is_axis == true)
	{
		int ne	= Solution.setup.NS+Solution.setup.ND;
#pragma omp parallel for num_threads(nt)
		for (int c = 0;  c <= grid_data.NCM-1; c++)
		{
			double	S	= grid_data.cells.data_ptr[c]->S;
			Solution.Rn[c][Solution.setup.NS+1]	-= S * Solution.Wc[c][ne] * fsgn(grid_data.cells.data_ptr[c]->x[1]);
		}
	}


	// Check error
#pragma omp parallel for num_threads(nt)
	for (int c = 0; c <= grid_data.NCM-1; c++)
	{
		for (int i = 0; i <= Solution.setup.VAR-1; i++)
		{
			if (error_check_double(Solution.Rn[c][i]))
			{
				Error_message_type	error_message;
				error_message.module_name	="CFD-Inviscid Flux Residue";
				error_message.location_primary_name	= "Cell";
				error_message.location_primary		= c;
				error_message.location_secondary_name	= "Variable";
				error_message.location_secondary		= i;
				error_message.message	= " Problem in inviscid flux residue at cell!";
				error_message.print_message();
			}
		}
	}
}

