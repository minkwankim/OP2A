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





void DerivativesType5::d2TdQ2(Data::DataStorage& data_Q, Data::DataStorage& data_V,  	Data::DataStorage& data_MIX,	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT, Data::DataStorage2D& d2T)
{
	double rho_Cvtr = data_MIX(4);
	double rho = data_MIX(0);


	double sum_u2	= 0.0;
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	sum_u2	+= pow(data_V(k), 2.0);
	sum_u2 /= rho;

	int index_e1 = species_set.NS+ND;
	int index_e2 = species_set.NS+ND;
	int NE	= 2;
	int VAR	= data_Q.numData;

	// A. Species, s1
#pragma ivdep
	for (int s1= 0; s1 <= species_set.NS-1; s1++)
	{
		// a. Species, s2
		for (int s2 = 0; s2 <= species_set.NS-1; s2++)
		{
			d2T(s1, s2)	= (-sum_u2 - species_set.species[s1].Cv_tra*dT(s2) - species_set.species[s2].Cv_tra*dT(s1)) / rho_Cvtr;
		}

		// b. Momentum, k2
		for (int k2 = species_set.NS; k2 <= species_set.NS+ND-1; k2++)
		{
			d2T(s1, k2)	= (data_V(k2)/rho - species_set.species[s1].Cv_tra*dT(k2)) / rho_Cvtr;
		}

		// c. Energy, s2
		for (int e2 = index_e2; e2 <= index_e2+NE-1; e2++)
		{
			d2T(s1, e2) = (-species_set.species[s1].Cv_tra*dT(e2)) / rho_Cvtr;
		}
	}

	for (int k1 = species_set.NS; k1 <= species_set.NS+ND-1; k1++)
	{
#pragma ivdep
		for (int s2 = 0; s2 <= species_set.NS-1; s2++)
		{
			d2T(k1, s2)	= (data_V(k1)/rho - species_set.species[s2].Cv_tra*dT(k1)) / rho_Cvtr;
		}

#pragma ivdep
		for (int k2 = species_set.NS; k2 <= species_set.NS+ND-1; k2++)
		{
			if (k1 == k2)	d2T(k1, k2)	= -(1.0/rho) / rho_Cvtr;
			else			d2T(k1, k2)	= 0.0;
		}

#pragma ivdep
		for (int e2 = index_e2; e2 <= index_e2+NE-1; e2++)
		{
			d2T(k1, e2) = 0.0;
		}
	}

	int e1	= index_e1;
#pragma ivdep
	for (int s2 = 0; s2 <= species_set.NS-1; s2++)
	{
		d2T(e1, s2)	= (-species_set.species[s2].Cv_tra*dT(e1)) / rho_Cvtr;
	}

#pragma ivdep
	for (int k2 = species_set.NS; k2 <= species_set.NS+ND-1; k2++)
	{
		d2T(e1, k2)	= 0.0;
	}

#pragma ivdep
	for (int e2 = index_e2; e2 <= index_e2+NE-1; e2++)
	{
		d2T(e1, e2) = 0.0;
	}

#pragma ivdep
	for (e1 = index_e1; e1 <= index_e1+NE-1; e1++)
	{
		for (int j = 0; j <= VAR-1; j++)
		{
			d2T(e1, j) = -d2T(index_e1, j);
		}
	}
}









}
}
