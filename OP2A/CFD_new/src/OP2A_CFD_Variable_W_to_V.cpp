/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Variable_W_to_V.cpp
 * 			-  
 *  
 */

#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"



void CFD_W_to_V(vector<double> W, vector<double> &V, CFD_variable_setup_ver2 setup, CFD_mixture_data variable_data, SPECIES_DATA_BASIC &species)
{
	V	= W;

	double p;
	double rho_R;
	if (setup.energy_flag == 1 || setup.energy_flag == 3 || setup.energy_flag == 5 || setup.energy_flag == 7)
	{
		double pe;
		double Te;

		Te = V[setup.NS+setup.ND + setup.ID_T[ELE]];
		pe = (V[setup.NS-1]*species.R[setup.NS-1]) * Te;

		p	= W[setup.NS+setup.ND]	- pe;
		rho_R	= variable_data.rho*variable_data.R_mix - (V[setup.NS-1]*species.R[setup.NS-1]);
	}
	else
	{
		p	= W[setup.NS+setup.ND];
		rho_R	= variable_data.rho*variable_data.R_mix;
	}


	double T = p / rho_R;
	V[setup.NS+setup.ND]	= T;
}
