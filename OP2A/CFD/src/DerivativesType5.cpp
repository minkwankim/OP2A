/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 26, 2015
 *      			Author: Minkwan Kim
 *
 * DerivativesType5.cpp
 * 			-  
 *  
 */


#include <omp.h>

#include "CFD/include/Derivatives.hpp"


namespace OP2A{
namespace CFD{



void DerivativesType5::dTdQ(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_MIX, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT, Data::DataStorage& dTr)
{
	double rho_Cvtra = data_MIX(4);
	double rho_Cvrot = data_MIX(5);

	double sum_u2	= 0.0;
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	sum_u2	+= pow(data_V(k), 2.0);
	sum_u2	*= 0.5;

	double T	= data_V(species_set.NS+ND);
	double Tr	= data_V(species_set.NS+ND+1);

#pragma ivdep
	for (int s= 0; s <= species_set.NS-1; s++)
	{
		dT(s)	= (-species_set.species[s].h0 + sum_u2 - species_set.species[s].Cv_tra*T) / rho_Cvtra;
		dTr(s)	= - species_set.species[s].Cv_rot*Tr / rho_Cvrot;
	}

#pragma ivdep
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		dT(k)	= -data_V(k) / rho_Cvtra;
		dTr(k)	= 0.0;
	}


	dT(species_set.NS+ND)	= 1.0 / rho_Cvtra;
	dTr(species_set.NS+ND)	= 0.0;

	dT(species_set.NS+ND+1)	 = -1.0 / rho_Cvtra;
	dTr(species_set.NS+ND+1) = 1.0 / rho_Cvrot;
}


void DerivativesType5::dpdQ(Data::DataStorage& data_V, Data::DataStorage& data_MIX,  Data::DataStorage& dT, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dp)
{
	double rho_R = data_MIX(1);
	double T	= data_V(species_set.NS+ND);

#pragma ivdep
	for (int s= 0; s <= species_set.NS-1; s++)
	{
		dp(s)	= species_set.species[s].R*T + rho_R*dT(s);
	}

#pragma ivdep
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		dp(k)	= rho_R * dT(k);
	}

	dp(species_set.NS+ND)	= rho_R * dT(species_set.NS+ND);
	dp(species_set.NS+ND+1)	= rho_R * dT(species_set.NS+ND+1);
}



double DerivativesType5::a2(Data::DataStorage& data_Q, Data::DataStorage& data_W, Data::DataStorage& data_MIX,  Data::DataStorage& dp, CHEM::SpeciesSet& species_set, int ND)
{
	double rho = data_MIX(0);

	double a2_rho = 0.0;

	for (int s= 0; s <= species_set.NS-1; s++)
	{
		a2_rho	+= data_Q(s)*dp(s);
	}

	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		a2_rho += data_Q(k)*dp(k);
	}

	a2_rho	+= (data_Q(species_set.NS+ND) + data_W(species_set.NS+ND)) * dp(species_set.NS+ND);
	a2_rho	+= data_Q(species_set.NS+ND+1) * dp(species_set.NS+ND+1);

	if (a2_rho < 0.0 || a2_rho == std::numeric_limits<double>::infinity() || a2_rho != a2_rho)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative value of a2: Need to check dp_dQ.");
	}

	//Common::ErrorCheckNonNegative<double>(a2_rho, "DerivativesType5: rho_a2 cannot be negative");
	return (a2_rho / rho);
}














}
}
