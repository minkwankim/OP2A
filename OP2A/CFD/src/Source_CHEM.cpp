/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 17, 2015
 *      			Author: Minkwan Kim
 *
 * Spurce_CHEM.cpp
 * 			-  
 *  
 */



#include <omp.h>
#include <time.h>

#include "CFD/include/Source_CHEM.hpp"
#include "Math/include/MathMisc.hpp"

namespace OP2A{
namespace CFD{

void Source_CHEM::S_chem(Data::DataStorageVector<Data::DataStorage>& data1D, CHEM::SpeciesSet& species_set, int ND,
				   	   	int type, unsigned int indexQ, unsigned int indexV, unsigned int indexW,
				   	   	Data::DataStorage& Sc)
{
	vector<double>	kf(species_set.NR, 0.0);
	vector<double>	kb(species_set.NR, 0.0);
	double Tt, Tr, Tv, Te;
	double Tcf, Tcb;

	int indexE	= species_set.NS + ND;

	switch(type)
	{
	case 1:
		Tt = data1D(indexV)(indexE);
		Tr = Tt;
		Tv = Tt;
		Te = Tt;
		break;

	case 2:
		Tt = data1D(indexV)(indexE);
		Tr = Tt;
		Tv = Tt;
		Te = data1D(indexV)(indexE+1);
		break;

	case 3:
		Tt = data1D(indexV)(indexE);
		Tr = Tt;
		Tv = data1D(indexV)(indexE+1);
		Te = Tt;
		break;

	case 4:
		Tt = data1D(indexV)(indexE);
		Tr = Tt;
		Tv = data1D(indexV)(indexE+1);
		Te = data1D(indexV)(indexE+2);
		break;

	case 5:
		Tt = data1D(indexV)(indexE);
		Tr = data1D(indexV)(indexE+1);
		Tv = Tt;
		Te = Tt;
		break;

	case 6:
		Tt = data1D(indexV)(indexE);
		Tr = data1D(indexV)(indexE+1);
		Tv = Tt;
		Te = data1D(indexV)(indexE+2);
		break;

	case 7:
		Tt = data1D(indexV)(indexE);
		Tr = data1D(indexV)(indexE+1);
		Tv = data1D(indexV)(indexE+2);
		Te = Tt;
		break;

	case 8:
		Tt = data1D(indexV)(indexE);
		Tr = data1D(indexV)(indexE+1);
		Tv = data1D(indexV)(indexE+2);
		Te = data1D(indexV)(indexE+3);
		break;
	}

	double n_mix = 0.0;
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		n_mix = data1D(indexQ)(s) / species_set.species[s].m;
	}

#pragma ivdep
	for (int k = 0; k <= species_set.NR-1; k++)
	{
		Tcf = CHEM::Calc_Tc(Tt, Tr, Tv, Te, static_cast<CHEM::ReactionType>(species_set.reactions[k].type), 0);
		Tcb = CHEM::Calc_Tc(Tt, Tr, Tv, Te, static_cast<CHEM::ReactionType>(species_set.reactions[k].type), 1);

		kf[k] = species_set.reactions[k].kf(Tcf);
		kb[k] = species_set.reactions[k].kb(Tcb, n_mix);
	}


	vector<double> rho_M(species_set.NS, 0.0);
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		rho_M[s] = data1D(indexQ)(s)/species_set.species[s].M * 1.0e-3;
	}


	vector<double> Rk(species_set.NR, 0.0);
	vector<double> Rk_f(species_set.NR, 1.0);
	vector<double> Rk_b(species_set.NR, 1.0);

#pragma ivdep
	for (int k = 0; k <= species_set.NR-1; k++)
	{
		for (int s = 0; s <= species_set.NS-1; s++)
		{
			if (species_set.reactions[k].alpha[s] != 0.0)	Rk_f[k] *= pow(rho_M[s], species_set.reactions[k].alpha[s]);
			if (species_set.reactions[k].beta[s]  != 0.0)	Rk_b[k] *= pow(rho_M[s], species_set.reactions[k].beta[s]);
		}

		Rk_f[k] *= kf[k];
		Rk_b[k] *= kb[k];
		Rk[k]	= Rk_f[k] - Rk_b[k];
	}


	vector<double> Ws(species_set.NS, 0.0);
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Ws[s] = 0.0;
		for (int k = 0; k <= species_set.NR-1; k++)
		{
			Ws[s]	+= species_set.reactions[k].beta_m_alpha[s]*Rk[k];
		}

		Sc(s) += Ws[s] * species_set.species[s].M;
	}
}




}
}
