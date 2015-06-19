/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 19, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChangeType5.cpp
 * 			-  
 *  
 */



#include <omp.h>

#include "CFD/include/VariableChangeTypes.hpp"


namespace OP2A{
namespace CFD{




void VariableChangeType5::V_to_Q(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;

	double T	= data_V(indexT);
	double Tr	= data_V(indexTr);

	double E, Er;


	// 1. Enthalpy of formation
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_V(s)*species_set.species[s].h0;


	// 2. Kinetic energy
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_V(s);

	double E_k = 0.0;
#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= rho_mix*data_V(k)*data_V(k);
	E_k *= 0.5;


	// 3. Translation
	double E_tra = 0.0;
#pragma omp parallel for reduction(+:E_tra)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_tra	+= data_V(s) * (species_set.species[s].Cv_tra*T);
	}


	// 4. Rotation
	Er = 0.0;
#pragma omp parallel for reduction(+:Er)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Er	+= data_V(s) * (species_set.species[s].Cv_rot*Tr);
	}


	E = E_tra + E_h0 + E_k + Er;

	data_Q(indexT)	= E;
	data_Q(indexTr)	= Er;
}


void VariableChangeType5::Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;

	double E	= data_Q(indexT);
	double Er	= data_Q(indexTr);

	double T, Tr;



	// Enthalpy of formation
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_Q(s)*species_set.species[s].h0;


	// Kinetic Energy
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);

	double E_k = 0.0;
#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_Q(k)*data_Q(k) / rho_mix;
	E_k *= 0.5;


	// 1. Cv_tra_mix
	double Cv_tra_mix = 0.0;
#pragma omp parallel for reduction(+:Cv_tra_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Cv_tra_mix	+= data_Q(s) * species_set.species[s].Cv_tra;
	}

	// 2. Cv_rot_mix
	double Cv_rot_mix = 0.0;
#pragma omp parallel for reduction(+:Cv_rot_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Cv_rot_mix	+= data_Q(s) * species_set.species[s].Cv_rot;
	}


	T	= (E - E_k - E_h0 - Er) / Cv_tra_mix;
	Tr	= Er / Cv_rot_mix;


	Common::ErrorCheckNonNegative<double>(T, "Temperature(Type5): Q => V");
	Common::ErrorCheckNonNegative<double>(Tr, "Rotational Temperature(Type5): Q => V");

	data_V(indexT)	= T;
	data_V(indexTr)	= Tr;
}




//////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * V ><=> W
 */

void VariableChangeType5::V_to_W(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;

	double T	= data_V(indexT);
	double Tr	= data_V(indexTr);

	double p;


	double rhoRmix	= 0.0;
#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_V(s)*species_set.species[s].R;

	p	= rhoRmix * T;

	// Pressure
	data_W(indexT)	= p;
	data_W(indexTr)	= Tr;
}

void VariableChangeType5::W_to_V(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;

	double p	= data_W(indexT);
	double Tr	= data_W(indexTr);
	double T;


	double rhoRmix	= 0.0;
#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_W(s)*species_set.species[s].R;
	T	= p / rhoRmix;

	data_V(indexT)	= T;
	data_V(indexTr)	= Tr;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Q <=> W
 */
void VariableChangeType5::Q_to_W(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;

	double E	= data_Q(indexT);
	double Er	= data_Q(indexTr);

	double T, Tr, p;



	// Enthalpy of formation
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_Q(s)*species_set.species[s].h0;


	// Kinetic Energy
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);

	double E_k = 0.0;
#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_Q(k)*data_Q(k) / rho_mix;
	E_k *= 0.5;


	// 1. Cv_tra_mix
	double Cv_tra_mix = 0.0;
#pragma omp parallel for reduction(+:Cv_tra_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Cv_tra_mix	+= data_Q(s) * species_set.species[s].Cv_tra;
	}

	// 2. Cv_rot_mix
	double Cv_rot_mix = 0.0;
#pragma omp parallel for reduction(+:Cv_rot_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Cv_rot_mix	+= data_Q(s) * species_set.species[s].Cv_rot;
	}


	T	= (E - E_k - E_h0 - Er) / Cv_tra_mix;
	Tr	= Er / Cv_rot_mix;

	Common::ErrorCheckNonNegative<double>(T, "Temperature(Type5): Q => V");
	Common::ErrorCheckNonNegative<double>(Tr, "Rotational Temperature(Type5): Q => V");

	double rhoRmix	= 0.0;
#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_Q(s)*species_set.species[s].R;
	p	= rhoRmix * T;

	data_W(indexT)	= p;
	data_W(indexTr)	= Tr;
}



void VariableChangeType5::W_to_Q(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;

	double p	= data_W(indexT);
	double Tr	= data_W(indexTr);
	double T;


	double rhoRmix	= 0.0;
#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_W(s)*species_set.species[s].R;
	T	= p / rhoRmix;



	double E, Er;


	// 1. Enthalpy of formation
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_W(s)*species_set.species[s].h0;


	// 2. Kinetic energy
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_W(s);

	double E_k = 0.0;
#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= rho_mix*data_W(k)*data_W(k);
	E_k *= 0.5;


	// 3. Translation
	double E_tra = 0.0;
#pragma omp parallel for reduction(+:E_tra)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_tra	+= data_W(s) * (species_set.species[s].Cv_tra*T);
	}


	// 4. Rotation
	Er = 0.0;
#pragma omp parallel for reduction(+:Er)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Er	+= data_W(s) * (species_set.species[s].Cv_rot*Tr);
	}


	E = E_tra + E_h0 + E_k + Er;


	data_Q(indexT)	= E;
	data_Q(indexTr)	= Er;

}




}
}
