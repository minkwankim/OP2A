/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 8, 2015
 *      			Author: Minkwan Kim
 *
 * Derivatives_da2dQ.cpp
 * 			-  
 *  
 */




#include <omp.h>

#include "Common/include/MultiDimension.hpp"
#include "CFD/include/Derivatives.hpp"


namespace OP2A{
namespace CFD{

void Derivatives::da2dQ(Data::DataStorageVector<Data::DataStorage>& data1D,
						Data::DataStorage& dp,
						Data::DataStorage2D& d2p,
						CHEM::SpeciesSet& species_set, int ND,
						unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX,
						Data::DataStorage& da2)
{
	int VAR	= data1D(indexQ).numData;
	int indexE = species_set.NS + ND;

	double rho = data1D(indexMIX)(0);

	vector <double>				Qp(VAR, 0.0);
	vector <vector<double> >	dQp	= Common::vector_2D(VAR, VAR, 0.0);

	// 1. Assign Qp
#pragma ivdep
	for (int i = 0; i <= VAR-1; i++)
	{
		Qp[i]	= data1D(indexQ)(i) / rho;
	}

	Qp[indexE]	+= data1D(indexW)(indexE) / rho;


	// 2. Calculae dQp
#pragma ivdep
	for (int s1 = 0; s1 <= species_set.NS-1; s1++)
	{
		dQp[s1][s1]	= 1.0;

		for (int s2 = 0; s2 <= species_set.NS-1; s2++)
		{
			dQp[s1][s2]	+= -Qp[s1];
		}

		for (int s2 = 0; s2 <= species_set.NS-1; s2++)
		{
			dQp[s1][s2]	/= rho;
		}
	}

	for (int k1 = species_set.NS; k1 <= species_set.NS+ND-1; k1++)
	{
		double rhouk_rho2	= Qp[k1] / rho;

#pragma ivdep
		for (int s2 = 0; s2 <= species_set.NS-1; s2++)
		{
			dQp[k1][s2]	= -rhouk_rho2;
		}

		dQp[k1][k1]	= 1.0/rho;
	}


#pragma ivdep
	for (int s2 = 0; s2 <= species_set.NS-1; s2++)
	{
		dQp[indexE][s2] = (dp(s2) - Qp[indexE]) / rho;
	}

#pragma ivdep
	for (int k2 = species_set.NS; k2 <= species_set.NS+ND-1; k2++)
	{
		dQp[indexE][k2]	= dp(k2) / rho;
	}

#pragma ivdep
	for (int e2 = indexE; e2 <= VAR-1; e2++)
	{
		dQp[indexE][e2]	= dp(e2) / rho;
	}

	dQp[indexE][indexE]	+= 1.0 / rho;

#pragma ivdep
	for (int e1 = indexE+1; e1 <= VAR-1; e1++)
	{
		double E_rho2	= Qp[e1] / rho;
		for (int s2 = 0; s2 <= species_set.NS-1; s2++)
		{
			dQp[e1][s2] = -E_rho2;
		}

		dQp[e1][e1]	= 1.0/rho;
	}


#pragma ivdep
	for (int j = 0; j <= VAR-1; j++)
	{
		double aux1 = 0.0;
		for (int i = 0; i <= VAR-1; i++)
		{
			aux1	+= dQp[i][j]*dp(i) + Qp[i]*d2p(i, j);
		}

		da2(j)	= aux1;
	}
}


void Derivatives::dadQ(Data::DataStorage& da2, double a, Data::DataStorage& da)
{
	int VAR = da2.numData;
	double two_a = 2.0 *a;

#pragma ivdep
	for (int i = 0; i <= VAR-1; i++)
	{
		da(i)	= da2(i) / two_a;
	}
}




}
}
