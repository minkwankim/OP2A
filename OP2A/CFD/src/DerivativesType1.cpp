/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 25, 2015
 *      			Author: Minkwan Kim
 *
 * DerivativesType1.cpp
 * 			-  
 *  
 */



#include <omp.h>

#include "CFD/include/Derivatives.hpp"


namespace OP2A{
namespace CFD{

void DerivativesType1::dTdQ(Data::DataStorage& data_Q, Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT)
{
	double rho_Cvtr = 0.0;

#pragma omp parallel for reduction(+:rho_Cvtr)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_Cvtr += data_Q(s)*species_set.species[s].Cv_tr;

	double sum_u2	= 0.0;
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	sum_u2	+= pow(data_V(k), 2.0);
	sum_u2	*= 0.5;

	double T	= data_V(species_set.NS+ND);
	for (int s= 0; s <= species_set.NS-1; s++)
		dT(s)	= (-species_set.species[s].h0 + sum_u2 - species_set.species[s].Cv_tr*T) / rho_Cvtr;

	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
		dT(k)	= data_V(k) / rho_Cvtr;

	dT(species_set.NS+ND)	= 1.0 / rho_Cvtr;
}

void DerivativesType1::dpdQ(Data::DataStorage& data_V, Data::DataStorage& dT, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dp)
{
	double rho_R = 0.0;

#pragma omp parallel for reduction(+:rho_R)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_R += data_V(s)*species_set.species[s].R;

	double T	= data_V(species_set.NS+ND);

	for (int s= 0; s <= species_set.NS-1; s++)
		dp(s)	= species_set.species[s].R*T + rho_R*dT(s);

	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
		dp(k)	= rho_R * dT(k);

	dp(species_set.NS+ND)	= rho_R * dT(species_set.NS+ND);
}


double DerivativesType1::a2(Data::DataStorage& data_Q, Data::DataStorage& data_W, Data::DataStorage& dp, CHEM::SpeciesSet& species_set, int ND)
{
	double rho = 0.0;

#pragma omp parallel for reduction(+:rho)
	for (int s = 0; s <= species_set.NS-1; s++)	rho += data_Q(s);


	double a2_rho = 0.0;

	for (int s= 0; s <= species_set.NS-1; s++)	a2_rho	+= data_Q(s)*dp(s);

	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	a2_rho += data_Q(k)*dp(k);

	a2_rho	+= (data_Q(species_set.NS+ND) + data_W(species_set.NS+ND)) * dp(species_set.NS+ND);


	Common::ErrorCheckNonNegative<double>(a2_rho, "DerivativesType1: rho_a2 cannot be negative");
	return (a2_rho / rho);
}




}
}
