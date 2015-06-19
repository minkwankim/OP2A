/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 18, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChangeCommon.cpp
 * 			-  
 *  
 */

#include <omp.h>

#include "CFD/include/VariableChangeCommon.hpp"


namespace OP2A{
namespace CFD{

//  Q and V
void VariableChangeCommon::V_to_Q(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	// Density
#pragma omp parallel for num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	data_Q(s)	= data_V(s);

	// Momentum
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_V(s);

#pragma omp parallel for num_threads(CFD_NT)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	data_Q(k)	= data_V(k)*rho_mix;
}

void VariableChangeCommon::Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	// Density
#pragma omp parallel for num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	data_V(s)	= data_Q(s);

	// Momentum
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);

#pragma omp parallel for num_threads(CFD_NT)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	data_V(k)	= data_Q(k)/rho_mix;
}



// V and W
void VariableChangeCommon::V_to_W(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	// Density
#pragma omp parallel for num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	data_W(s)	= data_V(s);

	// Momentum
#pragma omp parallel for num_threads(CFD_NT)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	data_W(k)	= data_V(k);
}

void VariableChangeCommon::W_to_V(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	// Density
#pragma omp parallel for num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	data_V(s)	= data_W(s);

	// Momentum
#pragma omp parallel for num_threads(CFD_NT)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	data_V(k)	= data_W(k);
}



//  Q and W
void VariableChangeCommon::W_to_Q(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	// Density
#pragma omp parallel for num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	data_Q(s)	= data_W(s);

	// Momentum
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_W(s);

#pragma omp parallel for num_threads(CFD_NT)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	data_Q(k)	= data_W(k)*rho_mix;
}

void VariableChangeCommon::Q_to_W(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	// Density
#pragma omp parallel for num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	data_W(s)	= data_Q(s);

	// Momentum
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);

#pragma omp parallel for num_threads(CFD_NT)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	data_W(k)	= data_Q(k)/rho_mix;
}




}
}
