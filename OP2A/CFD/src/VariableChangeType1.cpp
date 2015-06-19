/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 16, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChange_Type1.cpp
 * 			-  
 *  
 */


#include <omp.h>

#include "CFD/include/VariableChangeTypes.hpp"


namespace OP2A{
namespace CFD{



void VariableChangeType1::V_to_Q(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_V(s);


	// Total Energy
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_V(s)*species_set.species[s].h0;

	double E_k = 0.0;
#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= rho_mix*data_V(k)*data_V(k);
	E_k *= 0.5;


	int indexT	= species_set.NS + ND;
	double E_int = 0.0;
#pragma omp parallel for reduction(+:E_int)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_int	+= data_V(s) * (species_set.species[s].Cv_tr*data_V(indexT));
		//E_int	+= data_V(s) * (species_set.species[s].Cv_tr*data_V(indexT) + species_set.species[s].e_VIB(data_V(indexT)) + species_set.species[s].e_ELE(data_V(indexT)));
	}

	data_Q(indexT)	= E_int + E_k + E_h0;
}



void VariableChangeType1::Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);


	// Translational Temperature
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_Q(s)*species_set.species[s].h0;

	double E_k = 0.0;
#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_Q(k)*data_Q(k) / rho_mix;
	E_k *= 0.5;


	int indexT	= species_set.NS + ND;

	double Cv_mix = 0.0;
#pragma omp parallel for reduction(+:Cv_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Cv_mix	+= data_Q(s) * species_set.species[s].Cv_tr;
	}

	double T	= (data_Q(indexT) - E_k - E_h0) / Cv_mix;
	Common::ErrorCheckNonNegative<double>(T, "Temperature(Type1)");
	data_V(indexT)	= T;
}



//////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * V ><=> W
 */

void VariableChangeType1::V_to_W(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	double rhoRmix	= 0.0;
#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_V(s)*species_set.species[s].R;


	// Pressure
	int indexT	= species_set.NS + ND;
	data_W(indexT)	= rhoRmix*data_V(indexT);
}


void VariableChangeType1::W_to_V(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	double rhoRmix	= 0.0;
#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_W(s)*species_set.species[s].R;


	// Temperature
	int indexT	= species_set.NS + ND;
	data_V(indexT)	= data_W(indexT) / rhoRmix;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Q <=> W
 */
void VariableChangeType1::Q_to_W(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);


	// Translational Temperature
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_Q(s)*species_set.species[s].h0;

	double E_k = 0.0;
#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_Q(k)*data_Q(k) / rho_mix;
	E_k *= 0.5;


	int indexT	= species_set.NS + ND;
	double Cv_mix = 0.0;
#pragma omp parallel for reduction(+:Cv_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Cv_mix	+= data_Q(s) * species_set.species[s].Cv_tr;
	}

	double T	= (data_Q(indexT) - E_k - E_h0) / Cv_mix;


	double rhoRmix	= 0.0;
#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_Q(s)*species_set.species[s].R;


	// Pressure
	data_W(indexT)	= rhoRmix*T;
}


void VariableChangeType1::W_to_Q(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	double rhoRmix	= 0.0;
#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_W(s)*species_set.species[s].R;

	// Temperature
	int indexT	= species_set.NS + ND;
	double T = data_W(indexT) / rhoRmix;

	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_W(s);


	// Total Energy
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_W(s)*species_set.species[s].h0;

	double E_k = 0.0;
#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= rho_mix*data_W(k)*data_W(k);
	E_k *= 0.5;

	double E_int = 0.0;
#pragma omp parallel for reduction(+:E_int)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_int	+= data_W(s) * (species_set.species[s].Cv_tr*T);
	}

	data_Q(indexT)	= E_int + E_k + E_h0;
}






}
}
