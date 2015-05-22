/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 15, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Variable_convert_all.cpp
 * 			-  
 *  
 */


#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"
#include "../include/OP2A_CFD_Variable_change.hpp"

void CFD_Calculate_all_required_variables_inviscid(SOL_CFD &Solution, int NCM, SPECIES_DATA_BASIC	&species, int NT, int mode)
{
	switch (mode)
	{
	case 0:
#pragma omp parallel for num_threads(NT)
		for (int c = 0; c <= NCM-1; c++)
		{
			Solution.mixture_data_c[c].calculate_data(Solution.setup.NS, Solution.setup.ID_T[ROT], Solution.Qc[c], species.Cv_t, species.Cv_r, species.R, species.M);
			CFD_Q_to_V(Solution.Qc[c], Solution.Vc[c], Solution.setup, Solution.mixture_data_c[c], species);
			CFD_V_to_W(Solution.Vc[c], Solution.Wc[c], Solution.setup, Solution.mixture_data_c[c], species);
		}
		break;

	case 1:
#pragma omp parallel for num_threads(NT)
		for (int c = 0; c <= NCM-1; c++)
		{
			Solution.mixture_data_c[c].calculate_data(Solution.setup.NS, Solution.setup.ID_T[ROT], Solution.Vc[c], species.Cv_t, species.Cv_r, species.R, species.M);
			CFD_V_to_Q(Solution.Vc[c], Solution.Wc[c], Solution.setup, Solution.mixture_data_c[c], species);
			CFD_V_to_W(Solution.Vc[c], Solution.Wc[c], Solution.setup, Solution.mixture_data_c[c], species);
		}
		break;

	case 2:
#pragma omp parallel for num_threads(NT)
		for (int c = 0; c <= NCM-1; c++)
		{
			Solution.mixture_data_c[c].calculate_data(Solution.setup.NS, Solution.setup.ID_T[ROT], Solution.Wc[c], species.Cv_t, species.Cv_r, species.R, species.M);
			CFD_W_to_V(Solution.Wc[c], Solution.Vc[c], Solution.setup, Solution.mixture_data_c[c], species);
			CFD_V_to_Q(Solution.Vc[c], Solution.Qc[c], Solution.setup, Solution.mixture_data_c[c], species);
		}
		break;
	}
}
