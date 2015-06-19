/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 18, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChange.cpp
 * 			-  
 *  
 */



#include <omp.h>

#include "CFD/include/VariableChange.hpp"


namespace OP2A{
namespace CFD{


void VariableChange::V_to_Q(unsigned int type, unsigned int CFD_NT, Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q)
{
	VariableChangeCommon::V_to_Q(data_V, species_set, ND, data_Q, CFD_NT);

	switch (type)
	{
	case 1:
		VariableChangeType1::V_to_Q(data_V, species_set, ND, data_Q, CFD_NT);
		break;

	case 2:
		VariableChangeType2::V_to_Q(data_V, species_set, ND, data_Q, CFD_NT);
		break;
	case 3:
		VariableChangeType3::V_to_Q(data_V, species_set, ND, data_Q, CFD_NT);
		break;
	case 4:
		break;
	case 5:
		break;
	case 6:
		break;
	case 7:
		break;
	case 8:
		break;
	}
}






void VariableChange::Q_to_V(unsigned int variabletype, unsigned int CFD_NT, Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V)
{
	VariableChangeCommon::Q_to_V(data_Q, species_set, ND, data_V, CFD_NT);

	switch (variabletype)
	{
	case 1:
		VariableChangeType1::Q_to_V(data_Q, species_set, ND, data_V, CFD_NT);
		break;

	case 2:
		VariableChangeType2::Q_to_V(data_Q, species_set, ND, data_V, CFD_NT);
		break;
	case 3:
		VariableChangeType3::Q_to_V(data_Q, species_set, ND, data_V, CFD_NT);
		break;
	case 4:
		break;
	case 5:
		break;
	case 6:
		break;
	case 7:
		break;
	case 8:
		break;
	}
}



void VariableChange::V_to_W(unsigned int variabletype, unsigned int CFD_NT, Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W)
{
	VariableChangeCommon::V_to_W(data_V, species_set, ND, data_W, CFD_NT);

	switch (variabletype)
	{
	case 1:
		VariableChangeType1::V_to_W(data_V, species_set, ND, data_W, CFD_NT);
		break;

	case 2:
		VariableChangeType2::V_to_W(data_V, species_set, ND, data_W, CFD_NT);
		break;
	case 3:
		break;
	case 4:
		break;
	case 5:
		break;
	case 6:
		break;
	case 7:
		break;
	case 8:
		break;
	}
}


void VariableChange::W_to_V(unsigned int variabletype, unsigned int CFD_NT, Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V)
{
	VariableChangeCommon::W_to_V(data_W, species_set, ND, data_V, CFD_NT);

	switch (variabletype)
	{
	case 1:
		VariableChangeType1::W_to_V(data_W, species_set, ND, data_V, CFD_NT);
		break;

	case 2:
		VariableChangeType2::W_to_V(data_W, species_set, ND, data_V, CFD_NT);
		break;
	case 3:
		break;
	case 4:
		break;
	case 5:
		break;
	case 6:
		break;
	case 7:
		break;
	case 8:
		break;
	}
}


void VariableChange::W_to_Q(unsigned int variabletype, unsigned int CFD_NT, Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q)
{
	VariableChangeCommon::W_to_Q(data_W, species_set, ND, data_Q, CFD_NT);

	switch (variabletype)
	{
	case 1:
		VariableChangeType1::W_to_Q(data_W, species_set, ND, data_Q, CFD_NT);
		break;

	case 2:
		VariableChangeType2::W_to_Q(data_W, species_set, ND, data_Q, CFD_NT);
		break;
	case 3:
		break;
	case 4:
		break;
	case 5:
		break;
	case 6:
		break;
	case 7:
		break;
	case 8:
		break;
	}
}

void VariableChange::Q_to_W(unsigned int variabletype, unsigned int CFD_NT, Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W)
{
	VariableChangeCommon::Q_to_W(data_Q, species_set, ND, data_W, CFD_NT);

	switch (variabletype)
	{
	case 1:
		VariableChangeType1::Q_to_W(data_Q, species_set, ND, data_W, CFD_NT);
		break;

	case 2:
		VariableChangeType2::Q_to_W(data_Q, species_set, ND, data_W, CFD_NT);
		break;
	case 3:
		break;
	case 4:
		break;
	case 5:
		break;
	case 6:
		break;
	case 7:
		break;
	case 8:
		break;
	}
}






}
}
