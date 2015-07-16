/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 15, 2015
 *      			Author: Minkwan Kim
 *
 * Source_NONEQ_S_he.cpp
 * 			-  
 *  
 */




#include <omp.h>
#include <time.h>

#include "CFD/include/Source_NONEQ.hpp"
#include "Math/include/MathMisc.hpp"

namespace OP2A{
namespace CFD{


double Source_NONEQ::S_he(Data::DataStorageVector<Data::DataStorage>& data1D, CHEM::SpeciesSet& species_set, int ND,
				   	   	   unsigned int indexTe, unsigned int indexQ, unsigned int indexV, unsigned int indexW)
{
	double S = 0.0;

	double T  = data1D(indexV)(species_set.NS+ND);
	double Te = data1D(indexV)(indexTe);

	double Vth;
	double aux1;
	double sum = 0.0;


	double ne_cgs = 0.0;

#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			ne_cgs = data1D(indexQ)(s) / species_set.species[s].m/1.0e6;
			break;
		}
	}

	double temp1 = 8.0*MATH_PI/27.0 * pow(ELECTRON_CHARGE_CGS_2/C_BOLTZMANN_CGS/Te,2.0);
	double temp2 = 9.0 * pow(C_BOLTZMANN_CGS, 3.0) * pow(Te, 3.0) / (4.0*MATH_PI*pow(ELECTRON_CHARGE_CGS, 6.0));

	double sigma_en	= 2.0e-19;
	double sigma_ei	= temp1 * 1.e-4 * log(temp2);


	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			if(species_set.species[s].charge == 0.0)
			{
				sum	+= data1D(indexQ)(s)*N_AVOGADRO_SI/species_set.species[s].M/species_set.species[s].M * sigma_en;
			}
			else
			{
				sum	+= data1D(indexQ)(s)*N_AVOGADRO_SI/species_set.species[s].M/species_set.species[s].M * sigma_ei;
			}
		}
		else
		{
			aux1 	= 3*Ru*data1D(indexQ)(s) * (T - Te);
			Vth 	= sqrt(8.0*species_set.species[s].R*Te / MATH_PI);
		}
	}

	S = aux1 * Vth * sum;
	return (S);
}




}
}
