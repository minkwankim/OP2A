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
		VariableChangeType4::V_to_Q(data_V, species_set, ND, data_Q, CFD_NT);
		break;

	case 5:
		VariableChangeType5::V_to_Q(data_V, species_set, ND, data_Q, CFD_NT);
		break;

	case 6:
		VariableChangeType6::V_to_Q(data_V, species_set, ND, data_Q, CFD_NT);
		break;

	case 7:
		VariableChangeType7::V_to_Q(data_V, species_set, ND, data_Q, CFD_NT);
		break;

	case 8:
		VariableChangeType8::V_to_Q(data_V, species_set, ND, data_Q, CFD_NT);
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
		VariableChangeType4::Q_to_V(data_Q, species_set, ND, data_V, CFD_NT);
		break;

	case 5:
		VariableChangeType5::Q_to_V(data_Q, species_set, ND, data_V, CFD_NT);
		break;

	case 6:
		VariableChangeType6::Q_to_V(data_Q, species_set, ND, data_V, CFD_NT);
		break;

	case 7:
		VariableChangeType7::Q_to_V(data_Q, species_set, ND, data_V, CFD_NT);
		break;

	case 8:
		VariableChangeType8::Q_to_V(data_Q, species_set, ND, data_V, CFD_NT);
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
		VariableChangeType3::V_to_W(data_V, species_set, ND, data_W, CFD_NT);
		break;

	case 4:
		VariableChangeType4::V_to_W(data_V, species_set, ND, data_W, CFD_NT);
		break;

	case 5:
		VariableChangeType5::V_to_W(data_V, species_set, ND, data_W, CFD_NT);
		break;

	case 6:
		VariableChangeType6::V_to_W(data_V, species_set, ND, data_W, CFD_NT);
		break;

	case 7:
		VariableChangeType7::V_to_W(data_V, species_set, ND, data_W, CFD_NT);
		break;

	case 8:
		VariableChangeType8::V_to_W(data_V, species_set, ND, data_W, CFD_NT);
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
		VariableChangeType3::W_to_V(data_W, species_set, ND, data_V, CFD_NT);
		break;

	case 4:
		VariableChangeType4::W_to_V(data_W, species_set, ND, data_V, CFD_NT);
		break;

	case 5:
		VariableChangeType5::W_to_V(data_W, species_set, ND, data_V, CFD_NT);
		break;

	case 6:
		VariableChangeType6::W_to_V(data_W, species_set, ND, data_V, CFD_NT);
		break;

	case 7:
		VariableChangeType7::W_to_V(data_W, species_set, ND, data_V, CFD_NT);
		break;

	case 8:
		VariableChangeType8::W_to_V(data_W, species_set, ND, data_V, CFD_NT);
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
		VariableChangeType3::W_to_Q(data_W, species_set, ND, data_Q, CFD_NT);
		break;

	case 4:
		VariableChangeType4::W_to_Q(data_W, species_set, ND, data_Q, CFD_NT);
		break;

	case 5:
		VariableChangeType5::W_to_Q(data_W, species_set, ND, data_Q, CFD_NT);
		break;

	case 6:
		VariableChangeType6::W_to_Q(data_W, species_set, ND, data_Q, CFD_NT);
		break;

	case 7:
		VariableChangeType7::W_to_Q(data_W, species_set, ND, data_Q, CFD_NT);
		break;

	case 8:
		VariableChangeType8::W_to_Q(data_W, species_set, ND, data_Q, CFD_NT);
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
		VariableChangeType3::Q_to_W(data_Q, species_set, ND, data_W, CFD_NT);
		break;

	case 4:
		VariableChangeType4::Q_to_W(data_Q, species_set, ND, data_W, CFD_NT);
		break;

	case 5:
		VariableChangeType5::Q_to_W(data_Q, species_set, ND, data_W, CFD_NT);
		break;

	case 6:
		VariableChangeType6::Q_to_W(data_Q, species_set, ND, data_W, CFD_NT);
		break;

	case 7:
		VariableChangeType7::Q_to_W(data_Q, species_set, ND, data_W, CFD_NT);
		break;

	case 8:
		VariableChangeType8::Q_to_W(data_Q, species_set, ND, data_W, CFD_NT);
		break;
	}
}





void VariableChange::From_Q(unsigned int variabletype, Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
{
	switch (variabletype)
	{
	case 1:
		VariableChangeType1::From_Q(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 2:
		VariableChangeType2::From_Q(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 3:
		VariableChangeType3::From_Q(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 4:
		VariableChangeType4::From_Q(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 5:
		VariableChangeType5::From_Q(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 6:
		VariableChangeType6::From_Q(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 7:
		VariableChangeType7::From_Q(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 8:
		VariableChangeType8::From_Q(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;
	}
}


void VariableChange::From_V(unsigned int variabletype, Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
{
	switch (variabletype)
	{
	case 1:
		VariableChangeType1::From_V(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 2:
		VariableChangeType2::From_V(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 3:
		VariableChangeType3::From_V(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 4:
		VariableChangeType4::From_V(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 5:
		VariableChangeType5::From_V(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 6:
		VariableChangeType6::From_V(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 7:
		VariableChangeType7::From_V(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 8:
		VariableChangeType8::From_V(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	}
}

void VariableChange::From_W(unsigned int variabletype, Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
{
	switch (variabletype)
	{
	case 1:
		VariableChangeType1::From_W(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 2:
		VariableChangeType2::From_W(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 3:
		VariableChangeType3::From_W(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 4:
		VariableChangeType4::From_W(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 5:
		VariableChangeType5::From_W(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 6:
		VariableChangeType6::From_W(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 7:
		VariableChangeType7::From_W(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	case 8:
		VariableChangeType8::From_W(data_Q, data_V, data_W, data_MIX, data_Xs, data_Ys, species_set, ND, CFD_NT);
		break;

	}
}




}
}
