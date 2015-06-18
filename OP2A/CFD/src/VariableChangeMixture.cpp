/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 16, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChangeMixture.cpp
 * 			-  
 *  
 */



#include <omp.h>

#include "CFD/include/VariableChangeMixture.hpp"


namespace OP2A{
namespace CFD{

// Mole fraction
void VariableChangeMixture::Xs(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Xs)
{
#pragma omp parallel for num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS; s++)	data_Xs(s)	= data_Q(s) / species_set.species[s].m;


	// Total number of Moles
	double n_mix	= 0.0;
#pragma omp parallel for reduction(+:n_mix)
	for (int s = 0; s <= species_set.NS; s++)	n_mix	+= data_Xs(s);

#pragma omp parallel for num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS; s++)	data_Xs(s)	/= n_mix;
}


// Mass fraction
void VariableChangeMixture::Ys(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Ys)
{
	// Total number of Moles
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS; s++)	rho_mix	+= data_Q(s);

#pragma omp parallel for num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS; s++)	data_Ys(s) = data_Q(s)/rho_mix;
}





}
}
