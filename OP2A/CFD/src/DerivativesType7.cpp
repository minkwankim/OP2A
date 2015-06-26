/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 26, 2015
 *      			Author: Minkwan Kim
 *
 * DerivativesType7.cpp
 * 			-  
 *  
 */



#include <omp.h>

#include "CFD/include/Derivatives.hpp"


namespace OP2A{
namespace CFD{

void DerivativesType7::dTdQ(Data::DataStorage& data_Q, Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT, Data::DataStorage& dTr, Data::DataStorage& dTv)
{
	double T	= data_V(species_set.NS+ND);
	double Tr	= data_V(species_set.NS+ND+1);
	double Tv	= data_V(species_set.NS+ND+2);



	double rho_Cvtra = 0.0;
	double rho_Cvrot = 0.0;
	double rho_CvVE  = 0.0;


#pragma omp parallel for reduction(+:rho_Cvtra)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		rho_Cvtra += data_Q(s)*species_set.species[s].Cv_tra;
	}

#pragma omp parallel for reduction(+:rho_Cvrot)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		rho_Cvrot += data_Q(s)*species_set.species[s].Cv_rot;
	}

#pragma omp parallel for reduction(+:rho_CvVE)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		rho_CvVE += data_Q(s)*species_set.species[s].Cv_VE(Tv);
	}


	double sum_u2	= 0.0;
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	sum_u2	+= pow(data_V(k), 2.0);
	sum_u2	*= 0.5;


	for (int s= 0; s <= species_set.NS-1; s++)
	{
		dT(s)	= (-species_set.species[s].h0 + sum_u2 - species_set.species[s].Cv_tra*T) / rho_Cvtra;
		dTr(s)	= - species_set.species[s].Cv_rot*Tr / rho_Cvrot;
		dTv(s)	= -species_set.species[s].e_VE(Tv) / rho_CvVE;
	}

	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		dT(k)	= data_V(k) / rho_Cvtra;
		dTr(k)	= 0.0;
		dTv(k)	= 0.0;
	}

	dT(species_set.NS+ND)	= 1.0 / rho_Cvtra;
	dTr(species_set.NS+ND)	= 0.0;
	dTv(species_set.NS+ND)	= 0.0;

	dT(species_set.NS+ND+1)		= -1.0 / rho_Cvtra;
	dTr(species_set.NS+ND+1)	= 1.0 / rho_Cvrot;
	dTv(species_set.NS+ND+1)	= 0.0;

	dT(species_set.NS+ND+2)		= -1.0 / rho_Cvtra;
	dTr(species_set.NS+ND+2)	= 0.0;
	dTv(species_set.NS+ND+2)	= 1.0 / rho_CvVE;
}



void DerivativesType7::dpdQ(Data::DataStorage& data_V, Data::DataStorage& dT, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dp)
{
	double rho_R = 0.0;

#pragma omp parallel for reduction(+:rho_R)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_R += data_V(s)*species_set.species[s].R;

	double T	= data_V(species_set.NS+ND);

	for (int s= 0; s <= species_set.NS-1; s++)
	{
		dp(s)	= species_set.species[s].R*T + rho_R*dT(s);
	}

	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		dp(k)	= rho_R * dT(k);
	}

	dp(species_set.NS+ND)	= rho_R * dT(species_set.NS+ND);
	dp(species_set.NS+ND+1)	= rho_R * dT(species_set.NS+ND+1);
	dp(species_set.NS+ND+2)	= rho_R * dT(species_set.NS+ND+2);
}


double DerivativesType7::a2(Data::DataStorage& data_Q, Data::DataStorage& data_W, Data::DataStorage& dp, CHEM::SpeciesSet& species_set, int ND)
{
	double rho = 0.0;

#pragma omp parallel for reduction(+:rho)
	for (int s = 0; s <= species_set.NS-1; s++)	rho += data_Q(s);


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
	a2_rho	+= data_Q(species_set.NS+ND+1)*dp(species_set.NS+ND+1);
	a2_rho	+= data_Q(species_set.NS+ND+2)*dp(species_set.NS+ND+2);

	Common::ErrorCheckNonNegative<double>(a2_rho, "DerivativesType7: rho_a2 cannot be negative");
	return (a2_rho / rho);
}





}
}
