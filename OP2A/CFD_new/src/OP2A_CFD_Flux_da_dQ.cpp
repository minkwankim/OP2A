/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 15, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Flux_da_dQ.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_CFD_Flux.hpp"


void CFD_Calculate_da2_dQ(CFD_variable_setup_ver2 setup, vector<double> &Q, vector<double> &V, vector<double> &W, CFD_mixture_data &variable_data, vector<double> &dp_dQ, vector < vector<double> > &d2p_dQ2,	vector<double> &da2_dQ, SPECIES_DATA_BASIC &species)
{
	vector <double>				Qp(setup.VAR, 0.0);
	vector <vector<double> >	dQp	= vector_2D(setup.VAR, setup.VAR, 0.0);


	for (int s = 0; s <= setup.VAR-1; s++)	Qp[s]	= Q[s] / variable_data.rho;
	Qp[setup.NS+setup.ND]	+= W[setup.NS+setup.ND] / variable_data.rho;


	for (int s = 0; s <= setup.NS-1; s++)
	{
		dQp[s][s]	= 1.0;
		for (int r = 0; r <= setup.NS-1; r++)	dQp[s][r]	+= -Qp[s];
		for (int r = 0; r <= setup.NS-1; r++)	dQp[s][r]	= dQp[s][r] / variable_data.rho;
	}

	for (int k = 0; k <= setup.ND-1; k++)
	{
		for (int r = 0; r <= setup.NS-1; r++)	dQp[setup.NS+k][r]			= -Qp[setup.NS+k] / variable_data.rho;
		dQp[setup.NS+k][setup.NS+k]	= 1.0/variable_data.rho;
	}

	for (int r = 0; r <= setup.NS-1; r++)	dQp[setup.NS+setup.ND][r]					= (dp_dQ[r] - Qp[setup.NS+setup.ND]) / variable_data.rho;
	for (int r = 0; r <= setup.ND-1; r++)	dQp[setup.NS+setup.ND][setup.NS+r]			= dp_dQ[setup.NS+r] / variable_data.rho;
	for (int r = 0; r <= setup.NE-1; r++)	dQp[setup.NS+setup.ND][setup.NS+setup.ND+r]	= dp_dQ[setup.NS+setup.ND+r] / variable_data.rho;
	dQp[setup.NS+setup.ND][setup.NS+setup.ND]	+= 1.0 / variable_data.rho;

	for (int k = 1; k <= setup.NE-1; k++)
	{
		for (int r = 0; r <= setup.NS-1; r++)	dQp[setup.NS+setup.ND+k][r]	= -Qp[setup.NS+k]/variable_data.rho;
		dQp[setup.NS+setup.ND+k][setup.NS+setup.ND+k]	= 1.0/variable_data.rho;
	}



	for (int r = 0; r <= setup.VAR-1; r++)
	{
		double aux1 = 0.0;
		for (int s = 0; s <= setup.VAR-1; s++)	aux1	+= dQp[s][r]*dp_dQ[s] + Qp[s]*d2p_dQ2[s][r];

		da2_dQ[r]	= aux1;
	}
}


void CFD_Calculate_da_dQ(CFD_variable_setup_ver2 setup, vector<double> &Q, vector<double> &V, vector<double> &W, CFD_mixture_data &variable_data, vector<double> &da2_dQ, double a, vector<double> &da_dQ, SPECIES_DATA_BASIC &species)
{
	for (int i = 0; i <= setup.VAR-1; i++)	da_dQ[i]	= da2_dQ[i] / (2.0*a);
}


double CFD_Calculate_a2(CFD_variable_setup_ver2 setup, vector<double> &Q, vector<double> &W, CFD_mixture_data &variable_data, vector<double> &dp_dQ, SPECIES_DATA_BASIC &species)
{
	double H	= (Q[setup.NS+setup.ND]	+ W[setup.NS+setup.ND]) / variable_data.rho;
	double	a2 = 0.0;

	for (int s = 0; s <= setup.NS-1; s++) a2 += Q[s]/variable_data.rho * dp_dQ[s];
	for (int k = 0; k <= setup.ND-1; k++) a2 += Q[setup.NS+k]/variable_data.rho	* dp_dQ[setup.NS+k];
	a2 += H * dp_dQ[setup.NS+setup.ND];
	for (int m = 1; m <= setup.NE-1; m++) a2 += Q[setup.NS+setup.ND+m]/variable_data.rho *dp_dQ[setup.NS+setup.ND+m];

	return (a2);
}
