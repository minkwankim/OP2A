/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 4, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Flux_residue_AUSM.cpp
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


void CFD_residue_AUSM(GRID_CLASS &grid_data, SOL_CFD &Solution, SPECIES_DATA_BASIC	&species, bool is_axis, int nt, double M_inf)
{
	// Step 2: Calculate residues
#pragma omp parallel for num_threads(nt)
	for (int c = 0; c <= grid_data.NCM-1; c++)
	{
		vector < vector <double> >	dT_dQ = vector_2D(Solution.setup.NE, Solution.setup.VAR, 0.0);
		CFD_Calculate_dT_dQ(Solution.setup, Solution.Qc[c], Solution.Vc[c], Solution.mixture_data_c[c], dT_dQ, species);

		vector<double>	dp_dQ(Solution.setup.VAR, 0.0);
		CFD_Calculate_dp_dQ(Solution.setup, Solution.Vc[c], Solution.mixture_data_c[c], dT_dQ, dp_dQ,	species);

		double a2		= CFD_Calculate_a2(Solution.setup, Solution.Qc[c], Solution.Wc[c], Solution.mixture_data_c[c], dp_dQ, species);

		double a;
		if (a2 >= 0.0)	a	= sqrt(a2);
		else			program_error_type1("Error negative a2 LEFT [AUSM RESIDUE]: Please check dp and dT calculation!!");

		double q2	= 0.0;
		for (int k = 0; k <= Solution.setup.ND-1; k++)	q2	+= Solution.Vc[c][Solution.setup.NS+k]*Solution.Vc[c][Solution.setup.NS+k];
		double M2	= q2/a2;

		double eps	= fmin(1.0, fmax(M2, M_inf*M_inf));
		double aux1 = (1.0 - eps)*(Solution.mixture_data_c[c].gamma - 1.0) / a2;

	}
}
