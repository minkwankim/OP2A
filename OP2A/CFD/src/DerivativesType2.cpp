/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 25, 2015
 *      			Author: Minkwan Kim
 *
 * DerivativesType2.cpp
 * 			-  
 *  
 */



#include <omp.h>

#include "CFD/include/Derivatives.hpp"


namespace OP2A{
namespace CFD{

void DerivativesType2::dTdQ(Data::DataStorage& data_Q, Data::DataStorage& data_V,		Data::DataStorage& data_MIX,	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT,  Data::DataStorage& dTe)
{
	double rho_Cvtr_wo_e = data_MIX(4);

	double sum_u2	= 0.0;
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	sum_u2	+= pow(data_V(k), 2.0);
	sum_u2	*= 0.5;


	double T	= data_V(species_set.NS+ND);
	double Te	= data_V(species_set.NS+ND+1);
	double o_rhoe_Cve = 0.0;

#pragma ivdep
	for (int s= 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			dT(s)	= (-species_set.species[s].h0 + sum_u2 - species_set.species[s].Cv_tr*T) / rho_Cvtr_wo_e;
			dTe(s)	= 0.0;
		}
		else
		{
			o_rhoe_Cve = data_Q(s)*species_set.species[s].Cv_tra;
			if (o_rhoe_Cve != 0.0)	o_rhoe_Cve = 1.0 / o_rhoe_Cve;

			dT(s)	= sum_u2 / rho_Cvtr_wo_e;
			dTe(s)	= (-species_set.species[s].h0 - species_set.species[s].Cv_tra*Te) * o_rhoe_Cve;
		}
	}

#pragma ivdep
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		dT(k)	= -data_V(k) / rho_Cvtr_wo_e;
		dTe(k)	= 0.0;
	}

	dT(species_set.NS+ND)	= 1.0 / rho_Cvtr_wo_e;
	dTe(species_set.NS+ND)	= 0.0;

	dT(species_set.NS+ND+1)		= -1.0 / rho_Cvtr_wo_e;
	dTe(species_set.NS+ND+1)	= o_rhoe_Cve;
}



void DerivativesType2::dpdQ(Data::DataStorage& data_V, Data::DataStorage& data_MIX,	Data::DataStorage& dT, 			Data::DataStorage& dTe, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dp)
{
	double rho_R_wo_e = data_MIX(1);
	double T		= data_V(species_set.NS+ND);
	double Te		= data_V(species_set.NS+ND+1);
	double rhoe_Re 	= 0.0;


#pragma ivdep
	for (int s= 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			dp(s)	= species_set.species[s].R*T + rho_R_wo_e*dT(s);
		}
		else
		{
			rhoe_Re = data_V(s)*species_set.species[s].R;

			dp(s)	= rho_R_wo_e*dT(s) + species_set.species[s].R*Te + rhoe_Re*dTe(s);
		}
	}


#pragma ivdep
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		dp(k)	= rho_R_wo_e * dT(k);
	}

	dp(species_set.NS+ND)	= rho_R_wo_e*dT(species_set.NS+ND);
	dp(species_set.NS+ND+1)	= rho_R_wo_e*dT(species_set.NS+ND+1) + rhoe_Re*dTe(species_set.NS+ND+1);
}


double DerivativesType2::a2(Data::DataStorage& data_Q, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& dp, CHEM::SpeciesSet& species_set, int ND)
{
	double rho = data_MIX(0);

	double a2_rho = 0.0;

	for (int s= 0; s <= species_set.NS-1; s++)	a2_rho	+= data_Q(s)*dp(s);

	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	a2_rho += data_Q(k)*dp(k);

	a2_rho	+= (data_Q(species_set.NS+ND) + data_W(species_set.NS+ND)) * dp(species_set.NS+ND);
	a2_rho	+= data_Q(species_set.NS+ND+1)*dp(species_set.NS+ND+1);

	if (a2_rho < 0.0 || a2_rho == std::numeric_limits<double>::infinity() || a2_rho != a2_rho)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative value of a2(Type2): Need to check dp_dQ.");
	}
	//Common::ErrorCheckNonNegative<double>(a2_rho, "DerivativesType2: rho_a2 cannot be negative");
	return (a2_rho / rho);
}


}
}
