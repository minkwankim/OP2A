/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 15, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Flux_dT_dQ.cpp
 * 			-  
 *  
 */



#include "../include/OP2A_CFD_Flux.hpp"


void CFD_Calculate_dT_dQ(CFD_variable_setup_ver2 setup, vector<double> &Q, vector<double> &V, CFD_mixture_data &variable_data, vector< vector<double> > &dT_dQ, SPECIES_DATA_BASIC &species)
{
	CFD_Calculate_dTtra_dQ(setup, Q, V, variable_data, dT_dQ[0], species);

	switch (setup.energy_flag)
	{
	case 1:
		CFD_Calculate_dTe_dQ(setup, Q, V, variable_data, dT_dQ[1], species);
		break;

	case 2:
		CFD_Calculate_dTvib_dQ(setup, Q, V, variable_data, dT_dQ[1], species);
		break;

	case 3:
		CFD_Calculate_dTvib_dQ(setup, Q, V, variable_data, dT_dQ[1], species);
		CFD_Calculate_dTe_dQ(setup, Q, V, variable_data, dT_dQ[2], species);
		break;

	case 4:
		CFD_Calculate_dTrot_dQ(setup, Q, V, variable_data, dT_dQ[1], species);
		break;

	case 5:
		CFD_Calculate_dTrot_dQ(setup, Q, V, variable_data, dT_dQ[1], species);
		CFD_Calculate_dTe_dQ(setup, Q, V, variable_data, dT_dQ[2], species);
		break;

	case 6:
		CFD_Calculate_dTrot_dQ(setup, Q, V, variable_data, dT_dQ[1], species);
		CFD_Calculate_dTvib_dQ(setup, Q, V, variable_data, dT_dQ[2], species);
		break;

	case 7:
		CFD_Calculate_dTrot_dQ(setup, Q, V, variable_data, dT_dQ[1], species);
		CFD_Calculate_dTvib_dQ(setup, Q, V, variable_data, dT_dQ[2], species);
		CFD_Calculate_dTe_dQ(setup, Q, V, variable_data, dT_dQ[3], species);
		break;
	}
}

void CFD_Calculate_dTtra_dQ(CFD_variable_setup_ver2 setup, vector<double> &Q, vector<double> &V, CFD_mixture_data &variable_data, vector<double> &dT_dQ, SPECIES_DATA_BASIC &species)
{
	double rho_Cv	= variable_data.rho	* variable_data.Cv_bar;
	double alp	= 0.0;

	for (int k = 0; k <= setup.ND-1; k++)	alp += pow(V[setup.NS+k], 2.0);


	// Species
	if (setup.ID_T[ROT] == 0)
	{
		for (int s = 0; s <= setup.NS-1; s++)
			dT_dQ[s]	= (-species.h0[s] + 0.5*alp - V[setup.NS+setup.ND]*(species.Cv_t[s]+species.Cv_r[s])) / rho_Cv;
	}
	else
	{
		for (int s = 0; s <= setup.NS-1; s++)
			dT_dQ[s]	= (-species.h0[s] + 0.5*alp - V[setup.NS+setup.ND]*species.Cv_t[s]) / rho_Cv;
	}

	// Momentum
	for (int k = 0; k <= setup.ND-1; k++)	dT_dQ[setup.NS+k]	= -V[setup.NS+k] / rho_Cv;


	// Energy modes
	dT_dQ[setup.NS+setup.ND]	= 1.0 / rho_Cv;
	for (int m = 1; m <= setup.NE-1; m++)	dT_dQ[setup.NS+setup.ND+m]	-= 1.0/ rho_Cv;
}



void CFD_Calculate_dTrot_dQ(CFD_variable_setup_ver2 setup, vector<double> &Q, vector<double> &V, CFD_mixture_data &variable_data, vector<double> &dT_dQ, SPECIES_DATA_BASIC &species)
{
	double rho_Cv	= variable_data.rho	* variable_data.Cv_bar_rot;

	// Species
	int ne	= setup.NS+setup.ND;
	for (int s = 0; s <= setup.NS-1; s++)	dT_dQ[s]	= - V[ne+setup.ID_T[ROT]]*species.Cv_r[s] / rho_Cv;

	// Momentum
	for (int k = 0; k <= setup.ND-1; k++)	dT_dQ[setup.NS+k]	= 0.0;


	// Energy modes
	for (int m = 0; m <= setup.NE-1; m++)	dT_dQ[ne+m]	= 0.0;
	dT_dQ[ne+setup.ID_T[ROT]]	= 1.0 / rho_Cv;
}


void CFD_Calculate_dTvib_dQ(CFD_variable_setup_ver2 setup, vector<double> &Q, vector<double> &V, CFD_mixture_data &variable_data, vector<double> &dT_dQ, SPECIES_DATA_BASIC &species)
{
	int ne	= setup.NS + setup.ND;

	double Tv 		= V[ne + setup.ID_T[VIB]];
	double rho_Cv	= 0.0;

	for (int s = 0; s <= setup.NS-1; s++)
	{
		int ss = species.whereis[s];
		rho_Cv	+= V[s] * (*species.data_entire)[ss].Cv_VEE(Tv);

		dT_dQ[s]	= -(*species.data_entire)[ss].e_VEE(Tv);
	}


	// Species
	for (int s = 0; s <= setup.NS-1; s++)	dT_dQ[s]	/= rho_Cv;

	// Momentum
	for (int k = 0; k <= setup.ND-1; k++)	dT_dQ[setup.NS+k]	= 0.0;


	// Energy modes
	for (int m = 0; m <= setup.NE-1; m++)	dT_dQ[ne+m]	= 0.0;
	dT_dQ[ne+setup.ID_T[VIB]]	= 1.0 / rho_Cv;
}



void CFD_Calculate_dTe_dQ(CFD_variable_setup_ver2 setup, vector<double> &Q, vector<double> &V, CFD_mixture_data &variable_data, vector<double> &dT_dQ, SPECIES_DATA_BASIC &species)
{
	int ne	= setup.NS + setup.ND;

	double Te 		= V[ne + setup.ID_T[ELE]];
	double rho_Cv	= 0.0;

	rho_Cv	= V[setup.NS-1] * species.Cv_t[setup.NS-1];
	for (int s = 0; s <= setup.VAR-1; s++)	dT_dQ[s] = 0.0;


	// Species
	dT_dQ[setup.NS-1]			= -1.0 / rho_Cv;
	dT_dQ[ne+setup.ID_T[ELE]]	= 1.0 / rho_Cv;
}
