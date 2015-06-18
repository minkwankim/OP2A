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


void VariableChange::V_to_Q(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q)
{
	VariableChangeCommon::V_to_Q(data_V, species_set, ND, data_Q);

	switch (variableType)
	{
	case 1:
		VariableChangeType1::V_to_Q(data_V, species_set, ND, data_Q);
		break;

	case 2:
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

void VariableChange::Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V)
{
	VariableChangeCommon::Q_to_V(data_Q, species_set, ND, data_V);

	switch (variableType)
	{
	case 1:
		VariableChangeType1::Q_to_V(data_Q, species_set, ND, data_V);
		break;

	case 2:
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



void VariableChange::V_to_W(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W)
{
	VariableChangeCommon::V_to_W(data_V, species_set, ND, data_W);

	switch (variableType)
	{
	case 1:
		VariableChangeType1::V_to_W(data_V, species_set, ND, data_W);
		break;

	case 2:
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

void VariableChange::W_to_V(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V)
{
	VariableChangeCommon::W_to_V(data_W, species_set, ND, data_V);

	switch (variableType)
	{
	case 1:
		VariableChangeType1::W_to_V(data_W, species_set, ND, data_V);
		break;

	case 2:
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


void VariableChange::W_to_Q(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q)
{
	VariableChangeCommon::W_to_Q(data_W, species_set, ND, data_Q);

	switch (variableType)
	{
	case 1:
		VariableChangeType1::W_to_Q(data_W, species_set, ND, data_Q);
		break;

	case 2:
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

void VariableChange::Q_to_W(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W)
{
	VariableChangeCommon::Q_to_W(data_Q, species_set, ND, data_W);

	switch (variableType)
	{
	case 1:
		VariableChangeType1::Q_to_W(data_Q, species_set, ND, data_W);
		break;

	case 2:
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
