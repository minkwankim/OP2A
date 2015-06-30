/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 25, 2015
 *      			Author: Minkwan Kim
 *
 * DerivativesType4.cpp
 * 			-  
 *  
 */



#include <omp.h>

#include "CFD/include/Derivatives.hpp"


namespace OP2A{
namespace CFD{

void DerivativesType4::dTdQ(Data::DataStorage& data_Q, Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT,  Data::DataStorage& dTv, Data::DataStorage& dTe)
{

	double T	= data_V(species_set.NS+ND);
	double Tv	= data_V(species_set.NS+ND+1);
	double Te	= data_V(species_set.NS+ND+2);

	double rho_Cvtr_wo_e = 0.0;
#pragma omp parallel for reduction(+:rho_Cvtr_wo_e)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rho_Cvtr_wo_e += data_Q(s)*species_set.species[s].Cv_tr;
		}
	}


	double rho_CvVE = 0.0;
#pragma omp parallel for reduction(+:rho_CvVE)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rho_CvVE += data_Q(s)*species_set.species[s].Cv_VE(Tv);
		}
	}

	double sum_u2	= 0.0;
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	sum_u2	+= pow(data_V(k), 2.0);
	sum_u2	*= 0.5;


	double o_rhoe_Cve = 0.0;

	for (int s= 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			dT(s)	= (-species_set.species[s].h0 + sum_u2 - species_set.species[s].Cv_tr*T) / rho_Cvtr_wo_e;
			dTv(s)	= -species_set.species[s].e_VE(Tv) / rho_CvVE;
			dTe(s)	= 0.0;
		}
		else
		{
			o_rhoe_Cve = data_Q(s)*species_set.species[s].Cv_tra;
			if (o_rhoe_Cve != 0.0)	o_rhoe_Cve = 1.0 / o_rhoe_Cve;

			dT(s)	= sum_u2  / rho_Cvtr_wo_e;
			dTv(s)	= 0.0;
			dTe(s)	= (-species_set.species[s].h0 - species_set.species[s].Cv_tra*Te) * o_rhoe_Cve;
		}
	}

	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		dT(k)	= -data_V(k) / rho_Cvtr_wo_e;
		dTv(k)	= 0.0;
		dTe(k)	= 0.0;
	}


	dT(species_set.NS+ND)	= 1.0 / rho_Cvtr_wo_e;
	dTv(species_set.NS+ND)	= 0.0;
	dTe(species_set.NS+ND)	= 0.0;


	dT(species_set.NS+ND+1)		= -1.0 / rho_Cvtr_wo_e;
	dTv(species_set.NS+ND+1)	= 1.0 / rho_CvVE;
	dTe(species_set.NS+ND+1)	= 0.0;


	dT(species_set.NS+ND+2)		= -1.0 / rho_Cvtr_wo_e;
	dTv(species_set.NS+ND+2)	= 0.0;
	dTe(species_set.NS+ND+2)	= o_rhoe_Cve;
}



void DerivativesType4::dpdQ(Data::DataStorage& data_V, Data::DataStorage& dT, Data::DataStorage& dTe, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dp)
{
	double rho_R_wo_e = 0.0;

#pragma omp parallel for reduction(+:rho_R_wo_e)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rho_R_wo_e += data_V(s)*species_set.species[s].R;
		}
	}



	double T	= data_V(species_set.NS+ND);
	double Te	= data_V(species_set.NS+ND+2);
	double rhoe_Re = 0.0;


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


	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		dp(k)	= rho_R_wo_e * dT(k);
	}

	dp(species_set.NS+ND)	= rho_R_wo_e*dT(species_set.NS+ND);
	dp(species_set.NS+ND+1)	= rho_R_wo_e*dT(species_set.NS+ND+1);
	dp(species_set.NS+ND+2)	= rho_R_wo_e*dT(species_set.NS+ND+2) + rhoe_Re*dTe(species_set.NS+ND+2);
}


double DerivativesType4::a2(Data::DataStorage& data_Q, Data::DataStorage& data_W, Data::DataStorage& dp, CHEM::SpeciesSet& species_set, int ND)
{
	double rho = 0.0;

#pragma omp parallel for reduction(+:rho)
	for (int s = 0; s <= species_set.NS-1; s++)	rho += data_Q(s);


	double a2_rho = 0.0;

	for (int s= 0; s <= species_set.NS-1; s++)	a2_rho	+= data_Q(s)*dp(s);

	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	a2_rho += data_Q(k)*dp(k);

	a2_rho	+= (data_Q(species_set.NS+ND) + data_W(species_set.NS+ND)) * dp(species_set.NS+ND);
	a2_rho	+= data_Q(species_set.NS+ND+1)*dp(species_set.NS+ND+1);
	a2_rho	+= data_Q(species_set.NS+ND+2)*dp(species_set.NS+ND+2);


	Common::ErrorCheckNonNegative<double>(a2_rho, "DerivativesType4: rho_a2 cannot be negative");
	return (a2_rho / rho);
}

}
}

