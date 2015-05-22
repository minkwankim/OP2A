/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 16, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Jacobian_source.cpp
 * 			-  
 *  
 */

#include "../include/OP2A_CFD_Jacobian.hpp"



void CFD_Jacobian_Source_terms(GRID_CLASS &grid, SOL_CFD &Solution,	SPECIES_DATA_BASIC &species, REACTION_DATA_ver2 &reactions,
								vector<vector<double> > &kf, vector<vector<double> > &Rf, vector<vector<double> > &kb, vector<vector<double> > &Rb,
		bool is_axis, bool is_viscous, bool is_chemical, int nt)
{
#pragma omp parallel for num_threads(nt)
	for (int c = 0; c <= grid.NCM-1; c++)
	{
		// Initizliae variables
		for (int i = 0; i <= Solution.setup.VAR-1; i++)
			for (int j = 0; j <= Solution.setup.VAR-1; j++)	Solution.Jacobian_source[c][i][j]	= 0.0;
	}



	if (is_chemical == true)
	{
		CFD_Jacobian_chemistry_single_fluid(grid, Solution, kf, Rf, kb, Rb, reactions, Solution.Jacobian_source, is_axis, nt);
	}

	if (is_axis == true)
	{
#pragma omp parallel for num_threads(nt)
		for (int c = 0; c <= grid.NCM-1; c++)
		{
			CFD_Jacobian_axis_inviscid(Solution.setup, Solution.dp_dQ[c], Solution.Jacobian_source[c], grid.cells.data_ptr[c]->S);
		}
	}
}



