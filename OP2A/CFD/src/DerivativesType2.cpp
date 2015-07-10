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



void DerivativesType2::d2TdQ2(Data::DataStorage& data_Q, Data::DataStorage& data_V,  	Data::DataStorage& data_MIX,	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT, Data::DataStorage& dTe, Data::DataStorage2D& d2T, Data::DataStorage2D& d2Te)
{
	// Basic variables
	double rho_Cvtr = data_MIX(4);
	double rho = data_MIX(0);

	double rhoe_Cve;
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoe_Cve	= species_set.species[s].Cv_tra * data_Q(s);
			break;
		}
	}


	double sum_u2	= 0.0;
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
		sum_u2	+= pow(data_V(k), 2.0);

	sum_u2 /= rho;

	int index_e1 = species_set.NS+ND;
	int index_e2 = species_set.NS+ND;
	int NE	= 2;
	int VAR	= data_Q.numData;

	// A. Species
#pragma ivdep
	for (int s1= 0; s1 <= species_set.NS-1; s1++)
	{
		if (species_set.species[s1].type != CHEM::SpeciesType::Electron)
		{
			// a. Species, s2
			for (int s2 = 0; s2 <= species_set.NS-1; s2++)
			{
				if (species_set.species[s2].type != CHEM::SpeciesType::Electron)
				{
					d2T(s1, s2)		= (-sum_u2 - species_set.species[s2].Cv_tr*dT(s1) - species_set.species[s1].Cv_tr*dT(s2)) / rho_Cvtr;
					d2Te(s1, s2)	= 0.0;
				}
				else
				{
					d2T(s1, s2)	= (-sum_u2 - species_set.species[s1].Cv_tr*dT(s2)) / rho_Cvtr;
					d2Te(s1, s2)	= 0.0;
				}
			}

			// b. Momentum, k2
			for (int k2 = species_set.NS; k2 <= index_e2-1; k2++)
			{
				d2T(s1, k2)		= (data_V(k2)/rho - species_set.species[s1].Cv_tr*dT(k2)) / rho_Cvtr;
				d2Te(s1, k2)	= 0.0;
			}

			// c. Energy: total and electron
			for (int e2 = index_e2; e2 <= index_e2+NE-1; e2++)
			{
				d2T(s1, e2)		= (-species_set.species[s1].Cv_tr*dT(e2)) 	/ rho_Cvtr;
				d2Te(s1, e2)	= 0.0;
			}
		}
		else
		{
			// a. Species, s2
			for (int s2 = 0; s2 <= species_set.NS-1; s2++)
			{
				if (species_set.species[s2].type != CHEM::SpeciesType::Electron)
				{
					d2T(s1, s2)		= (-sum_u2 - species_set.species[s2].Cv_tr*dT(s1)) / rho_Cvtr;
					d2Te(s1, s2) 	= 0.0;
				}
				else
				{
					d2T(s1, s2)		= (-sum_u2) / rho_Cvtr;
					d2Te(s1, s2)	= (-species_set.species[s2].Cv_tra*dTe(s1) - species_set.species[s1].Cv_tra*dTe(s2)) / rhoe_Cve;
				}
			}

			// b. Momentum, k2
			for (int k2 = species_set.NS; k2 <= index_e2-1; k2++)
			{
				d2T(s1, k2)		= (data_V(k2)/rho) / rho_Cvtr;
				d2Te(s1, k2)	= 0.0;
			}

			// c. Energy: total and electron
			for (int e2 = index_e2; e2 <= index_e2+NE-1; e2++)
			{
				d2T(s1, e2) 	= 0.0;
				d2Te(s1, e2)	= (-species_set.species[s1].Cv_tra * dTe(e2)) / rhoe_Cve;
			}
		}
	}


	// B. Momentum, k1
	for (int k1 = species_set.NS; k1 <= index_e1-1; k1++)
	{
		double uk_rho = data_V(k1)/rho;

		// a. Species, s2
#pragma ivdep
		for (int s2 = 0; s2 <= species_set.NS-1; s2++)
		{
			if (species_set.species[s2].type != CHEM::SpeciesType::Electron)
			{
				d2T(k1, s2)		= (uk_rho - species_set.species[s2].Cv_tr*dT(k1)) / rho_Cvtr;
				d2Te(k1, s2)	= 0.0;
			}
			else
			{
				d2T(k1, s2)		= (uk_rho) / rho_Cvtr;
				d2Te(k1, s2)	= 0.0;
			}
		}

		// b. Momentum, k2
#pragma ivdep
		for (int k2 = species_set.NS; k2 <= index_e2-1; k2++)
		{
			if (k1 == k2)	d2T(k1, k2)	= -(1.0/rho) / rho_Cvtr;
			else			d2T(k1, k2)	= 0.0;

			d2Te(k1, k2)	= 0.0;
		}


		// c. Energy: total and electron
#pragma ivdep
		for (int e2 = index_e2; e2 <= index_e2+NE-1; e2++)
		{
			d2T(k1, e2) 	= 0.0;
			d2Te(k1, e2)	= 0.0;
		}
	}


	// C. Energy
	//		C1. Total energy
	// a. Species, s2
#pragma ivdep
	for (int s2 = 0; s2 <= species_set.NS-1; s2++)
	{
		if (species_set.species[s2].type != CHEM::SpeciesType::Electron)
		{
			d2T(index_e1, s2)	= (-species_set.species[s2].Cv_tr*dT(index_e1)) / rho_Cvtr;
			d2Te(index_e1, s2)	= 0.0;
		}
		else
		{
			d2T(index_e1, s2)	= 0.0;
			d2Te(index_e1, s2)	= 0.0;
		}
	}

	// b. Momentum, k2
#pragma ivdep
	for (int k2 = species_set.NS; k2 <= index_e2-1; k2++)
	{
		d2T(index_e1, k2)	= 0.0;
		d2Te(index_e1, k2)	= 0.0;
	}

	// c. Energy. e2
#pragma ivdep
	for (int e2 = index_e2; e2 <= index_e2+NE-1; e2++)
	{
		d2T(index_e1, e2)	= 0.0;
		d2Te(index_e1, e2)	= 0.0;
	}

	//		C2. OTHER energy
#pragma ivdep
	for (int e1 = 1; e1 <= index_e1+NE-2; e1++)
	{
		for (int j = 0; j <= VAR-1; j++)
		{
			d2T(e1, j)	= -d2T(index_e1, j);
			d2Te(e1, j)	= 0.0;
		}
	}

	//		C3. Electron energy
	// a. Species, s2
	int e1e = index_e1+NE-1;
#pragma ivdep
	for (int s2 = 0; s2 <= species_set.NS-1; s2++)
	{
		if (species_set.species[s2].type != CHEM::SpeciesType::Electron)
		{
			d2T(e1e, s2)	= -d2T(index_e1, s2);
			d2Te(e1e, s2)	= 0.0;
		}
		else
		{
			d2T(e1e, s2)	= 0.0;
			d2Te(e1e, s2)	= (-species_set.species[s2].Cv_tra*dTe(e1e)) / rhoe_Cve;
		}
	}

	// b. Momentum, k2
#pragma ivdep
	for (int k2 = species_set.NS; k2 <= index_e2-1; k2++)
	{
		d2T(e1e, k2)	= 0.0;
		d2Te(e1e, k2)	= 0.0;
	}

	// c. Energy
#pragma ivdep
	for (int e2 = index_e2; e2 <= index_e2+NE-1; e2++)
	{
		d2T(e1e, e2)	= 0.0;
		d2Te(e1e, e2)	= 0.0;
	}
}



}
}
