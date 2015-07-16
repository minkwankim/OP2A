/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 15, 2015
 *      			Author: Minkwan Kim
 *
 * JacobianSourceTermS_he.cpp
 * 			-  
 *  
 */



#include <omp.h>
#include <time.h>

#include "CFD/include/Jacobians.hpp"
#include "Math/include/OP2A_Matrix.hpp"
#include "Math/include/MathMisc.hpp"
#include "Common/include/Time_Info.hpp"
#include "Common/include/MultiDimension.hpp"

namespace OP2A{
namespace CFD{

void Jacobians::S_he(Data::DataStorageVector<Data::DataStorage>& data1D, Data::DataStorage2D& dT, CHEM::SpeciesSet& species_set, int ND,
					unsigned int indexTe, unsigned int indexQ, unsigned int indexV, unsigned int indexW,
					double S, Data::DataStorage& dS)
{
	int index_dTe = indexTe - species_set.NS - ND;
	double T  = data1D(indexV)(species_set.NS+ND);
	double Te = data1D(indexV)(indexTe);


	double Ru_rhoe = 0.0;
	double sqrt_Vth = 0.0;
	double ne_cgs = 0.0;
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			ne_cgs = data1D(indexQ)(s) / species_set.species[s].m/1.0e6;

			Ru_rhoe = Ru * data1D(indexQ)(s);
			sqrt_Vth = sqrt(8.0*Ru/MATH_PI/species_set.species[s].M);
			break;
		}
	}

	double temp1 = 8.0*MATH_PI/27.0 * pow(ELECTRON_CHARGE_CGS_2/C_BOLTZMANN_CGS/Te,2.0);
	double temp2 = 9.0 * pow(C_BOLTZMANN_CGS, 3.0) * pow(Te, 3.0) / (4.0*MATH_PI*pow(ELECTRON_CHARGE_CGS, 6.0));

	double sigma_en	= 2.0e-19;
	double sigma_ei	= temp1 * 1.e-4 * log(temp2);


	double A;
	double B;
	double C = 0.0;

	int VAR = data1D(indexQ).numData;
	vector<double> dA_dQ (VAR, 0.0);
	vector<double> dB_dQ (VAR, 0.0);
	vector<double> dC_dQ (VAR, 0.0);


	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			if(species_set.species[s].charge == 0.0)
			{
				C	+= data1D(indexQ)(s)*N_AVOGADRO_SI/species_set.species[s].M/species_set.species[s].M * sigma_en;
				dC_dQ[s] = N_AVOGADRO_SI/species_set.species[s].M/species_set.species[s].M * sigma_en;
			}
			else
			{
				C	+= data1D(indexQ)(s)*N_AVOGADRO_SI/species_set.species[s].M/species_set.species[s].M * sigma_ei;
				dC_dQ[s] = N_AVOGADRO_SI/species_set.species[s].M/species_set.species[s].M * sigma_ei;
			}
		}
		else
		{
			A 	= 3*Ru*data1D(indexQ)(s) * (T - Te);
			B 	= sqrt(8.0*species_set.species[s].R*Te / MATH_PI);

			dA_dQ[s] += 3.0*Ru*(T - Te);
		}
	}



	double aux1 = 3.0 * Ru_rhoe;
	double aux2 = sqrt_Vth / (2.0*sqrt(Te));

#pragma ivdep
	for (int i = 0; i <= VAR-1; i++)
	{
		dA_dQ[i] += aux1 * (dT(0, i) - dT(index_dTe, i));
		dB_dQ[i] = aux2 * dT(0, i);
	}


#pragma ivdep
	for (int i = 0; i <= VAR-1; i++)
	{
		dS(i)	= (dA_dQ[i]*B*C + dB_dQ[i]*A*C + dC_dQ[i]*A*B) * S;
	}
}





}
}
