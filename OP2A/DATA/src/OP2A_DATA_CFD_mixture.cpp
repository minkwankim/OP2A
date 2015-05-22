/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_DATA_CFD_mixture.hpp
 * 			-  
 *  
 */

#include "../include/OP2A_DATA_CFD_mixture.hpp"
#include "../../CFD_new/include/OP2A_CFD_utilities.hpp"
#include "../../PHYSICS/include/constant_text.hpp"
#include "../../PHYSICS/include/mixture_properties.hpp"


CFD_mixture_data::CFD_mixture_data()
{
	rho			= 0.0;
	Cv_bar		= 0.0;
	Cv_bar_rot	= 0.0;
	R_mix		= 0.0;
	M_mix		= 0.0;
	gamma		= 0.0;
}

CFD_mixture_data::~CFD_mixture_data()
{

}

void CFD_mixture_data::calculate_data(int NS, int NER, vector<double> Q, vector<double> Cv_tra_s, vector<double> Cv_rot_s, vector<double> R_s, vector<double> M_s)
{
	double M_tot	= 0.0;
	Xs.resize(NS);
	Ys.resize(NS);


	rho	= 0.0;
	for (int s = 0; s <= NS-1; s++)	rho	+= Q[s];

	for (int s = 0; s <= NS-1; s++)
	{
		Xs[s]	= Q[s]/M_s[s];
		M_tot	+= Xs[s];
	}

	if (M_tot != 0.0)
	{
#pragma omp parallel for
		for (int i = 0; i <= NS-1; i++)
		{
			Xs[i]	/= M_tot;
			Ys[i]	= Q[i]	/ rho;
		}
	}
	else
	{
#pragma omp parallel for
		for (int i = 0; i <= NS-1; i++)
		{
			Xs[i]	= 0.0;
			Ys[i]	= 0.0;
		}
	}


	R_mix		= cal_mixture_properties(R_s, NS, Ys, Xs, GAS_CONSTANT);
	M_mix		= cal_mixture_properties(M_s, NS, Ys, Xs, MOLECULAR_WEIGHT);
	Cv_bar		= cal_mixture_properties(Cv_tra_s, NS, Ys, Xs, SPECIFIC_HEAT);
	Cv_bar_rot	= cal_mixture_properties(Cv_rot_s, NS, Ys, Xs, SPECIFIC_HEAT);

	if (NER	== 0)
	{
		Cv_bar		+= Cv_bar_rot;
		Cv_bar_rot	= 0.0;
	}

	if (Cv_bar != 0.0)	gamma	= 1.0 + R_mix / (Cv_bar + Cv_bar_rot);
	else				gamma	= 5.0/3.0;
}

