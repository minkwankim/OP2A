/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Variable_V_to_Q.cpp
 * 			-  
 *  
 */

#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"



void CFD_V_to_Q(vector<double> V, vector<double> &Q, CFD_variable_setup_ver2 setup, CFD_mixture_data variable_data, SPECIES_DATA_BASIC &species)
{
	double T;
	double rho_Cv;


	for (int s = 0; s <= setup.NS-1; s++)	Q[s]			= V[s];
	for (int k = 0; k <= setup.ND-1; k++)	Q[setup.NS+k]	= V[setup.NS+k] * variable_data.rho;


	int ne	= setup.NS + setup.ND;

	switch(setup.energy_flag)
	{
	case 1:
		Q[ne+1]	= (species.Cv_t[setup.NS-1]*Q[setup.NS-1]) * V[ne+1];
		break;

	case 2:
		Q[ne+1]	= 0.0;
		for (int s = 0; s <= setup.NS-1; s++)
		{
			int ss = species.whereis[s];
			Q[ne+1]	+=	Q[s] * (*species.data_entire)[ss].e_VEE(V[ne+1]);
		}
		break;

	case 3:
		Q[ne+1]	= 0.0;
		for (int s = 0; s <= setup.NS-1; s++)
		{
			int ss = species.whereis[s];
			Q[ne+1]	+=	Q[s] * (*species.data_entire)[ss].e_VEE(V[ne+1]);
		}

		Q[ne+2]	= (species.Cv_t[setup.NS-1]*Q[setup.NS-1]) * V[ne+2];
		break;

	case 4:
		Q[ne+1]	= (variable_data.rho*variable_data.Cv_bar_rot) * V[ne+1];
		break;

	case 5:
		Q[ne+1]	= (variable_data.rho*variable_data.Cv_bar_rot) * V[ne+1];
		Q[ne+2]	= (species.Cv_t[setup.NS-1]*Q[setup.NS-1]) * V[ne+2];
		break;

	case 6:
		Q[ne+1]	= (variable_data.rho*variable_data.Cv_bar_rot) * V[ne+1];
		Q[ne+2]	= 0.0;
		for (int s = 0; s <= setup.NS-1; s++)
		{
			int ss = species.whereis[s];
			Q[ne+2]	+=	Q[s] * (*species.data_entire)[ss].e_VEE(V[ne+2]);
		}
		break;

	case 7:
		Q[ne+1]	= (variable_data.rho*variable_data.Cv_bar_rot) * V[ne+1];
		Q[ne+2]	= 0.0;
		for (int s = 0; s <= setup.NS-1; s++)
		{
			int ss = species.whereis[s];
			Q[ne+2]	+=	Q[s] * (*species.data_entire)[ss].e_VEE(V[ne+2]);
		}
		Q[ne+3]	= (species.Cv_t[setup.NS-1]*Q[setup.NS-1]) * V[ne+3];
		break;
	}

	rho_Cv	= variable_data.rho * variable_data.Cv_bar;
	Q[ne]	= rho_Cv * V[ne];
	for (int s = 0; s <= setup.NS-1; s++)	Q[ne]	+= Q[s]*species.h0[s];
	for (int k = 0; k <= setup.ND-1; k++)	Q[ne]	+= 0.5*variable_data.rho*V[setup.NS+k]*V[setup.NS+k];
	for (int k = 1; k <= setup.NE-1; k++)	Q[ne] 	+= Q[setup.NS+setup.ND+k];
}
