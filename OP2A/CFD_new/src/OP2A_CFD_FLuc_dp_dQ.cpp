/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 15, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_FLuc_dp_dQ.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_CFD_Flux.hpp"


void CFD_Calculate_dp_dQ_type1(CFD_variable_setup_ver2 setup, vector<double> &V, CFD_mixture_data &variable_data, vector<vector<double> > &dT_dQ, vector<double> &dp_dQ,	SPECIES_DATA_BASIC &species)
{
	double rho_R	= variable_data.rho * variable_data.R_mix;

	for (int i = 0; i <= setup.VAR-1; i++)	dp_dQ[i]	= rho_R * dT_dQ[0][i];
	for (int s = 0; s <= setup.NS-1; s++)	dp_dQ[s]	+= species.R[s]*V[setup.NS+setup.ND];
}

void CFD_Calculate_dp_dQ_type2(CFD_variable_setup_ver2 setup, vector<double> &V, CFD_mixture_data &variable_data, vector<vector<double> > &dT_dQ, vector<double> &dp_dQ,	SPECIES_DATA_BASIC &species)
{
	double rho_R	= 0.0;
	double rhoe_Re	= 0.0;

	for (int s = 0; s <= setup.NS-1; s++)
	{
		int s_ID	= species.whereis[s];
		if ((*species.data_entire)[s_ID].basic_data.type == ELECTRON)	rhoe_Re	= V[s]*species.R[s];
		else															rho_R	+= V[s]*species.R[s];
	}


	for (int i = 0; i <= setup.VAR-1; i++)	dp_dQ[i]	= rho_R*dT_dQ[0][i]	+ rhoe_Re*dT_dQ[setup.ID_T[ELE]][i];

	int ne	= setup.NS+setup.ND;
	for (int s = 0; s <= setup.NS-1; s++)
	{
		int s_ID	= species.whereis[s];
		if ((*species.data_entire)[s_ID].basic_data.type == ELECTRON)	dp_dQ[s]	+= species.R[s]*V[ne];
		else															dp_dQ[s]	+= species.R[s]*V[ne+setup.ID_T[ELE]];
	}
}



void CFD_Calculate_dp_dQ(CFD_variable_setup_ver2 setup, vector<double> &V, CFD_mixture_data &variable_data, vector<vector<double> > &dT_dQ, vector<double> &dp_dQ,	SPECIES_DATA_BASIC &species)
{

	if (setup.energy_flag == 1 || setup.energy_flag == 3 || setup.energy_flag == 5 || setup.energy_flag == 7)
	{
		CFD_Calculate_dp_dQ_type2(setup, V, variable_data, dT_dQ, dp_dQ,	species);
	}
	else
	{
		CFD_Calculate_dp_dQ_type1(setup, V, variable_data, dT_dQ, dp_dQ,	species);
	}
}



void CFD_Calculate_d2p_dQ2(CFD_variable_setup_ver2 setup, vector<double> &V, CFD_mixture_data &variable_data, vector < vector<double> > &d2p_dQ2,	SPECIES_DATA_BASIC &species)
{
	double R_bar_Cv_bar	= variable_data.R_mix / variable_data.Cv_bar;

	double u_square	= 0.0;
	for (int k = 0; k <= setup.ND-1; k++)	u_square	+= V[setup.NS+k]*V[setup.NS+k];


	// Variable 1
	vector<double>	dR_Cv_dQ(setup.VAR, 0.0);
	double rho_Cv_bar2	= variable_data.rho * variable_data.Cv_bar * variable_data.Cv_bar;

	if (setup.ID_T[ROT]	== 0)	for (int s = 0; s <= setup.NS-1; s++) dR_Cv_dQ[s]	= (species.R[s]*variable_data.Cv_bar - variable_data.R_mix*(species.Cv_r[s]+species.Cv_t[s])) / rho_Cv_bar2;
	else						for (int s = 0; s <= setup.NS-1; s++) dR_Cv_dQ[s]	= (species.R[s]*variable_data.Cv_bar - variable_data.R_mix*(species.Cv_t[s])) / rho_Cv_bar2;



	// Variable 2
	vector<double>	u_square_dQ(setup.VAR, 0.0);
	double temp1	= -2.0*u_square/variable_data.rho;
	for (int s = 0; s <= setup.NS-1; s++)	u_square_dQ[s]			= temp1;
	for (int k = 0; k <= setup.ND-1; k++)	u_square_dQ[setup.NS+k]	= 2.0*V[setup.NS+k]/variable_data.rho;

	vector <vector<double> >	du_dQ	= vector_2D(setup.ND, setup.VAR, 0.0);
	for (int k = 0; k <= setup.ND-1; k++)
	{
		temp1	= V[setup.NS+k]/variable_data.rho;

		for (int s = 0; s <= setup.NS-1; s++)	du_dQ[k][s]				= -temp1;
		du_dQ[k][setup.ND+k]	= 1.0/variable_data.rho;
	}


	/*
	 * Calculate d2p_dQ2
	 */
	for (int s = 0; s <= setup.NS-1; s++)
	{
		temp1	= 0.5*u_square - species.h0[s];
		for (int k = 0; k <= setup.VAR-1; k++)
		{
			d2p_dQ2[s][k]	= dR_Cv_dQ[k]*temp1 + 0.5*R_bar_Cv_bar*u_square_dQ[k];
		}
	}

	for (int k = 0; k <= setup.ND-1; k++)
	{
		for (int i = 0; i <= setup.VAR-1; i++)
		{
			d2p_dQ2[setup.NS+k][i]	= -(dR_Cv_dQ[i]*V[setup.NS+k] + R_bar_Cv_bar*du_dQ[k][i]);
		}
	}

	for (int i = 0; i <= setup.VAR-1; i++)	d2p_dQ2[setup.NS+setup.ND][i]	= dR_Cv_dQ[i];
	for (int k = 1; k <= setup.NE-1; k++)
	{
		for (int i = 0; i <= setup.VAR-1; i++)	d2p_dQ2[setup.NS+setup.ND+k][i]	= -dR_Cv_dQ[i];
	}
}
