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

void DerivativesType1::dTdQ(Data::DataStorage& data_Q, Data::DataStorage& data_V,  Data::DataStorage& data_MIX, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT)
{
	double rho_Cvtr = data_MIX(4);

	double sum_u2	= 0.0;
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	sum_u2	+= pow(data_V(k), 2.0);
	sum_u2	*= 0.5;

	double T	= data_V(species_set.NS+ND);

#pragma ivdep
	for (int s= 0; s <= species_set.NS-1; s++)
	{
		dT.data[s]	= (-species_set.species[s].h0 + sum_u2 - species_set.species[s].Cv_tr*T) / rho_Cvtr;
	}

#pragma ivdep
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		dT.data[k]	= -data_V.data[k] / rho_Cvtr;
	}

	dT.data[species_set.NS+ND]	= 1.0 / rho_Cvtr;
}

void DerivativesType1::dpdQ(Data::DataStorage& data_V, Data::DataStorage& data_MIX,	Data::DataStorage& dT, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dp)
{
	double rho_R = data_MIX(1);
	double T	= data_V(species_set.NS+ND);

#pragma ivdep
	for (int s= 0; s <= species_set.NS-1; s++)
	{
		dp.data[s]	= species_set.species[s].R*T + rho_R*dT.data[s];
	}

#pragma ivdep
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		dp.data[k]	= rho_R * dT.data[k];
	}

	dp.data[species_set.NS+ND]	= rho_R * dT.data[species_set.NS+ND];
}


double DerivativesType1::a2(Data::DataStorage& data_Q, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& dp, CHEM::SpeciesSet& species_set, int ND)
{
	double rho = data_MIX(0);
	double a2_rho = 0.0;

//#pragma omp parallel for reduction(+:a2_rho) num_threads(2)
	/*
	 * Do not vectorize for small number of sum. It will be slower than before
	 */
	for (int s= 0; s <= species_set.NS-1; s++)	a2_rho	+= data_Q.data[s]*dp.data[s];

	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	a2_rho += data_Q.data[k]*dp.data[k];

	a2_rho	+= (data_Q.data[species_set.NS+ND] + data_W.data[species_set.NS+ND]) * dp.data[species_set.NS+ND];

	if (a2_rho < 0.0 || a2_rho == std::numeric_limits<double>::infinity() || a2_rho != a2_rho)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative value of a2: Need to check dp_dQ.");
	}
	//Common::ErrorCheckNonNegative<double>(a2_rho, "DerivativesType1: rho_a2 cannot be negative");
	return (a2_rho / rho);
}


void DerivativesType1::d2TdQ2(Data::DataStorage& data_Q, Data::DataStorage& data_V,  	Data::DataStorage& data_MIX,	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT, Data::DataStorage2D& dT2)
{
	double rho_Cvtr = data_MIX(4);
	double rho = data_MIX(0);


	double sum_u2	= 0.0;
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	sum_u2	+= pow(data_V(k), 2.0);
	sum_u2 /= rho;

	int e1 = species_set.NS+ND;
	int e2 = species_set.NS+ND;

	for (int s1= 0; s1 <= species_set.NS-1; s1++)
	{
#pragma ivdep
		for (int s2 = 0; s2 <= species_set.NS-1; s2++)
		{
			dT2(s1, s2)	= (-sum_u2 - species_set.species[s1].Cv_tr*dT(s2) - species_set.species[s2].Cv_tr*dT(s1)) / rho_Cvtr;
		}

#pragma ivdep
		for (int k2 = species_set.NS; k2 <= species_set.NS+ND-1; k2++)
		{
			dT2(s1, k2)	= (data_V(k2)/rho - species_set.species[s1].Cv_tr*dT(k2)) / rho_Cvtr;
		}

		dT2(s1, e2) = (-species_set.species[s1].Cv_tr*dT(e2)) / rho_Cvtr;
	}

	for (int k1 = species_set.NS; k1 <= species_set.NS+ND-1; k1++)
	{
#pragma ivdep
		for (int s2 = 0; s2 <= species_set.NS-1; s2++)
		{
			dT2(k1, s2)	= (data_V(k1)/rho - species_set.species[s2].Cv_tr*dT(k1)) / rho_Cvtr;
		}

#pragma ivdep
		for (int k2 = species_set.NS; k2 <= species_set.NS+ND-1; k2++)
		{
			if (k1 == k2)	dT2(k1, k2)	= -(1.0/rho) / rho_Cvtr;
			else			dT2(k1, k2)	= 0.0;
		}

		dT2(k1, e2) = 0.0;
	}


#pragma ivdep
	for (int s2 = 0; s2 <= species_set.NS-1; s2++)
	{
		dT2(e1, s2)	= (-species_set.species[s2].Cv_tr*dT(e1)) / rho_Cvtr;
	}

#pragma ivdep
	for (int k2 = species_set.NS; k2 <= species_set.NS+ND-1; k2++)
	{
		dT2(e1, k2)	= 0.0;
	}

	dT2(e1, e2) = 0.0;
}




}
}
