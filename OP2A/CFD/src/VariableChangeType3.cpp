/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 19, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChangeType3.cpp
 * 			-  
 *  
 */



#include <omp.h>

#include "CFD/include/VariableChangeTypes.hpp"


namespace OP2A{
namespace CFD{



void VariableChangeType3::V_to_Q(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTve	= indexT + 1;

	double T	= data_V(indexT);
	double Tve	= data_V(indexTve);


	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_h0	+= data_V(s)*species_set.species[s].h0;
	}


	// 2. E_k (Kinetic energy)
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_V(s);

	double E_k = 0.0;
#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		E_k	+= rho_mix*data_V(k)*data_V(k);
	}
	E_k *= 0.5;


	// 3. E_int (translation/rotation energy) [Except electrons]
	double E_int = 0.0;
#pragma omp parallel for reduction(+:E_int)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_int	+= data_V(s) * (species_set.species[s].Cv_tr*T);
	}


	// 4. E_VE (VIBRATION-ELECTRONIC energy)
	double E_ve = 0.0;
#pragma omp parallel for reduction(+:E_int)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_ve	+= data_V(s) * species_set.species[s].e_VE(Tve);
	}


	data_Q(indexT)	= E_int + E_k + E_h0 + E_ve;
	data_Q(indexTve)	= E_ve;
}




void VariableChangeType3::Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);


	int indexT		= species_set.NS + ND;
	int indexTve	= indexT + 1;

	double E	= data_Q(indexT);
	double Eve  = data_Q(indexTve);

	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_h0	+= data_Q(s)*species_set.species[s].h0;
	}


	// 2. E_k (Kinetic energy)
	double E_k = 0.0;
#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		E_k	+= data_Q(k)*data_Q(k) / rho_mix;
	}
	E_k *= 0.5;


	// 3. rhoCvmix
	double rhoCvmix = 0.0;
#pragma omp parallel for reduction(+:rhoCvmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		rhoCvmix	+= data_Q(s) * species_set.species[s].Cv_tr;
	}


	// 4. rhoCve
	double rhoCve = 0.0;
#pragma omp parallel for reduction(+:rhoCve)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		rhoCve	+= data_Q(s) * species_set.species[s].Cv_tra;
	}

	double T, Tve;
	T	= (E - E_h0 - E_k - Eve) / rhoCvmix;
	Tve	= CFD_calculate_Tve(data_V, species_set, Eve, 10000, 0.0001, CFD_NT);

	Common::ErrorCheckNonNegative<double>(T,   "Temperature(Type3): Q => V");
	Common::ErrorCheckNonNegative<double>(Tve, "Vibrational Temperature(Type3): Q => V");

	data_V(indexT)		= T;
	data_V(indexTve)	= Tve;
}



}
}
