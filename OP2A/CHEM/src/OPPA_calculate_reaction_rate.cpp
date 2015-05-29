/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 7, 2015
 *      			Author: Minkwan Kim
 *
 * OPPA_calculate_reaction_rate.cpp
 * 			-  
 *  
 */




#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <omp.h>

#include "../include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"


void calculate_reaction_rate(REACTION_DATA_ver2 &reaction_data, vector<double> &rhos, vector<double> &T, vector<double> &kf, vector<double> &kb, vector<double> &Rf, vector<double> &Rb)
{
	kmp_set_defaults("KMP_AFFINITY = scatter");


	// 1. Initialize Rf and Rb
	Rf.reserve(reaction_data.NR);
	Rb.reserve(reaction_data.NR);


	// 2. Calculate 10^-3 * rho_s / Ms
	vector<double>	rhos_Ms(reaction_data.NS, 0.0);

#pragma omp parallel for
	for (int s = 0; s <= reaction_data.NS-1; s++)	rhos_Ms[s]	= 0.001 * rhos[s]/reaction_data.species_data[s].basic_data.M;

	double n_mix	= 0.0;
	for (int s = 0; s <= reaction_data.NS-1; s++)	n_mix	= n_mix + rhos[s]/reaction_data.species_data[s].basic_data.m;
	n_mix	= n_mix * 1.0e-6;



	//3. Calculate Forward/backward reaction rate
	for (int k = 0; k <= reaction_data.NR-1; k++)
	{
		double Tf, Tb;
		calculate_reaction_temperature(reaction_data.reaction_k[k], T, Tf, Tb);
		kf[k]	= cal_kf(reaction_data.reaction_k[k], Tf);
		kb[k]	= cal_kb(reaction_data.reaction_k[k], n_mix, Tb);


		double Rf_temp	= 1.0;
		for (int j = 0; j <= reaction_data.NS-1; j++)	Rf_temp	*= pow(rhos_Ms[j], reaction_data.reaction_k[k].Reactant_coeff[j]);
		Rf[k]	= 1000.0*kf[k] * Rf_temp;

		double Rb_temp	= 1.0;
		for (int j = 0; j <= reaction_data.NS-1; j++)	Rf_temp	*= pow(rhos_Ms[j], reaction_data.reaction_k[k].Product_coeff[j]);
		Rb[k]	= 1000.0*kb[k] * Rb_temp;
	}
}




void calculate_production_rate(REACTION_DATA_ver2 &reaction_data, vector<double> &Rf, vector<double> &Rb, vector<double> &S_chem)
{
	kmp_set_defaults("KMP_AFFINITY = scatter");

	if (S_chem.size() != reaction_data.NS)	S_chem.resize(reaction_data.NS);

#pragma omp parallel for
	for (int s = 0; s <= reaction_data.NS-1; s++)	S_chem[s]	= 0.0;


#pragma omp parallel for
	for (int s = 0; s <= reaction_data.NS-1; s++)
	{
		for (int k = 0; k <= reaction_data.NR-1; k++)
		{
			double Rf_m_Rb	= Rf[k] - Rb[k];
			S_chem[s]	+= (reaction_data.reaction_k[k].Product_coeff[s] - reaction_data.reaction_k[k].Reactant_coeff[s]) * Rf_m_Rb;
		}
	}


#pragma omp parallel for
	for (int s = 0; s <= reaction_data.NS-1; s++)
		S_chem[s]	= S_chem[s] * reaction_data.species_data[s].basic_data.M;
}
