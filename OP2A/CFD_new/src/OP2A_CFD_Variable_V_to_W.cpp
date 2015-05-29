/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Variable_V_toW.cpp
 * 			-  
 *  
 */

#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"



void CFD_V_to_W(vector<double> V, vector<double> &W, CFD_variable_setup_ver2 setup, CFD_mixture_data variable_data, SPECIES_DATA_BASIC &species)
{
	W	= V;

	double p = 0.0;
	double T = 0.0;
	double Te = 0.0;

	T	= V[setup.NS+setup.ND];
	p	= variable_data.rho * variable_data.R_mix * T;

	if (setup.energy_flag == 1 || setup.energy_flag == 3 || setup.energy_flag == 5 || setup.energy_flag == 7)
	{
		Te = V[setup.NS+setup.ND + setup.ID_T[ELE]];
		p = p + (V[setup.NS-1]*species.R[setup.NS-1]) * (Te - T);
	}

	W[setup.NS+setup.ND]	= p;
}
