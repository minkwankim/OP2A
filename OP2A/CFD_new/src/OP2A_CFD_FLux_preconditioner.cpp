/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 5, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_FLux_preconditioner.cpp
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


void CFD_Flux_Preconditioner(GRID_CLASS &grid_data,
							SOL_CFD &Solution,
							SPECIES_DATA_BASIC	&species,
							bool is_axis, int nt, double M_inf)
{
	int ne	= Solution.setup.NS + Solution.setup.ND;

	// Step 2: Calculate residues
#pragma omp parallel for num_threads(nt)
	for (int c = 0; c <= grid_data.NCM-1; c++)
	{
		CFD_Calculate_dT_dQ(Solution.setup, Solution.Qc[c], Solution.Vc[c], Solution.mixture_data_c[c], Solution.dT_dQ[c], species);
		CFD_Calculate_dp_dQ(Solution.setup, Solution.Vc[c], Solution.mixture_data_c[c], Solution.dT_dQ[c], Solution.dp_dQ[c], species);

		double a2	= CFD_Calculate_a2(Solution.setup, Solution.Qc[c], Solution.Wc[c], Solution.mixture_data_c[c], Solution.dp_dQ[c], species);

		double a;
		if (a2 >= 0.0)	a	= sqrt(a2);
		else			program_error_type1("Error negative a2 LEFT [AUSM RESIDUE]: Please check dp and dT calculation!!");

		double q2	= 0.0;
		for (int k = 0; k <= Solution.setup.ND-1; k++)	q2	+= Solution.Vc[c][Solution.setup.NS+k]*Solution.Vc[c][Solution.setup.NS+k];
		double M2	= q2/a2;

		double H	= (Solution.Qc[c][ne] + Solution.Wc[c][ne]) / Solution.mixture_data_c[c].rho;

		double eps	= fmin(1.0, fmax(M2, M_inf*M_inf));
		double aux1 = -(1.0 - eps)*(Solution.mixture_data_c[c].gamma - 1.0) / a2;

		for (int i = 0; i <= Solution.setup.VAR-1; i++)	Solution.Preconditioner[c][i][i] = 1.0;

		for (int s = 0; s <= Solution.setup.NS-1; s++)
		{
			//double Ys	= Solution.Qc[s] / Solution.mixture_data_c[c].rho;
			//for (int s1 = 0; s1 <= Solution.setup.NS-1; s1++)	Solution.Preconditioner[c][s][s1]					+= aux1*Ys * q2/2.0;
			//for (int k  = 0; s  <= Solution.setup.ND-1; k++)	Solution.Preconditioner[c][s][Solution.setup.NS+k]	+= aux1*Ys * -(Solution.Vc[c][Solution.setup.NS+k]);
			//Solution.Preconditioner[c][s][ne]	=
			//for (int k  = 0; s  <= Solution.setup.ND-1; k++)	Solution.Preconditioner[c][s][Solution.setup.NS+k]	+= aux1*Ys * -(Solution.Vc[c][Solution.setup.NS+k]);


		}

	}
}
