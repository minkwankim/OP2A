/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 22, 2015
 *      			Author: Minkwan Kim
 *
 * VariableCellData1D.cpp
 * 			-  
 *  
 */



#include "CFD/include/Variables.hpp"


namespace OP2A{
namespace CFD{

void CellData1D::CompleteDataWIGCase1(Data::DataStorageVector<Data::DataStorage>& cell1D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, int indexQ, int indexV, int indexW, int indexMIX, int indexXs, int indexYs)
{
	// 1. Calculate V
	CFD::VariableChange::Q_to_V(variabletype, CFD_NT, cell1D(indexQ), species_set, ND, cell1D(indexV));

	// 2. Calculate W
	CFD::VariableChange::V_to_W(variabletype, CFD_NT, cell1D(indexV), species_set, ND, cell1D(indexW));

	// 3. Calculate MIX
	CFD::VariableChangeMixture::MIX(cell1D(indexQ), species_set, ND, variabletype, CFD_NT, cell1D(indexYs), cell1D(indexXs), cell1D(indexMIX));
}


void CellData1D::CompleteDataWIGCase2(Data::DataStorageVector<Data::DataStorage>& cell1D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, int indexQ, int indexV, int indexW, int indexMIX, int indexXs, int indexYs)
{
	// 1. Calculate Q
	CFD::VariableChange::V_to_Q(variabletype, CFD_NT, cell1D(indexV), species_set, ND, cell1D(indexQ));

	// 2. Calculate W
	CFD::VariableChange::V_to_W(variabletype, CFD_NT, cell1D(indexV), species_set, ND, cell1D(indexW));

	// 3. Calculate MIX
	CFD::VariableChangeMixture::MIX(cell1D(indexQ), species_set, ND, variabletype, CFD_NT, cell1D(indexYs), cell1D(indexXs), cell1D(indexMIX));
}


void CellData1D::CompleteDataWIGCase3(Data::DataStorageVector<Data::DataStorage>& cell1D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, int indexQ, int indexV, int indexW, int indexMIX, int indexXs, int indexYs)
{
	// 1. Calculate V
	CFD::VariableChange::W_to_V(variabletype, CFD_NT, cell1D(indexW), species_set, ND, cell1D(indexV));

	// 2. Calculate Q
	CFD::VariableChange::V_to_Q(variabletype, CFD_NT, cell1D(indexV), species_set, ND, cell1D(indexQ));

	// 3. Calculate MIX
	CFD::VariableChangeMixture::MIX(cell1D(indexQ), species_set, ND, variabletype, CFD_NT, cell1D(indexYs), cell1D(indexXs), cell1D(indexMIX));
}


void CellData1D::CompleteDataWIGCase4(Data::DataStorageVector<Data::DataStorage>& cell1D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, int indexQ, int indexV, int indexW, int indexMIX, int indexXs, int indexYs)
{
	// 1. Calculate W
	CFD::VariableChange::V_to_W(variabletype, CFD_NT, cell1D(indexV), species_set, ND, cell1D(indexW));

	// 2. Calculate MIX
	CFD::VariableChangeMixture::MIX(cell1D(indexQ), species_set, ND, variabletype, CFD_NT, cell1D(indexYs), cell1D(indexXs), cell1D(indexMIX));
}



void CellData1D::InitializeOtherData(Data::DataStorageVector<Data::DataStorage>& cell1D, unsigned int CFD_NT, int indexR, int indexdQ, int indexQnew, int indexS, int indexDivV, int indexdpdQ, int indexDs,int indexMus)
{
	int numVar;

	// 1. Residue
	numVar	= cell1D(indexR).numData;
#pragma omp parallel for num_threads(CFD_NT)
	for (int i = 0; i <= numVar-1; i++)	cell1D(indexR).data[i]	= 0.0;

	// 2. dQ
	numVar	= cell1D(indexdQ).numData;
#pragma omp parallel for num_threads(CFD_NT)
	for (int i = 0; i <= numVar-1; i++)	cell1D(indexdQ).data[i]	= 0.0;

	// 3. Qnew
	numVar	= cell1D(indexQnew).numData;
#pragma omp parallel for num_threads(CFD_NT)
	for (int i = 0; i <= numVar-1; i++)	cell1D(indexQnew).data[i]	= 0.0;

	// 4. Source
	numVar	= cell1D(indexS).numData;
#pragma omp parallel for num_threads(CFD_NT)
	for (int i = 0; i <= numVar-1; i++)	cell1D(indexS).data[i]	= 0.0;

	// 5. DivV
	if (indexDivV != -1)
	{
		numVar	= cell1D(indexDivV).numData;
#pragma omp parallel for num_threads(CFD_NT)
		for (int i = 0; i <= numVar-1; i++)	cell1D(indexDivV).data[i]	= 0.0;
	}

	// 6. dpdQ
	if (indexdpdQ != -1)
	{
		numVar	= cell1D(indexdpdQ).numData;
#pragma omp parallel for num_threads(CFD_NT)
		for (int i = 0; i <= numVar-1; i++)	cell1D(indexdpdQ).data[i]	= 0.0;
	}


	// 7. Ds
	if (indexDs != -1)
	{
		numVar	= cell1D(indexDs).numData;
#pragma omp parallel for num_threads(CFD_NT)
		for (int i = 0; i <= numVar-1; i++)	cell1D(indexDs).data[i]	= 0.0;
	}

	// 8 Mus
	if (indexMus != -1)
	{
		numVar	= cell1D(indexMus).numData;
#pragma omp parallel for num_threads(CFD_NT)
		for (int i = 0; i <= numVar-1; i++)	cell1D(indexMus).data[i]	= 0.0;
	}
}




void CellData1D::CompleteDataWIG(Data::DataStorageVector<Data::DataStorage>& cell1D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT,
							int indexQ, int indexV, int indexW, int indexMIX, int indexXs, int indexYs,
							int indexR, int indexdQ, int indexQnew, int indexS, int indexDivV, int indexdpdQ, int indexDs,int indexMus,
							int typeCase, bool is_initialize)
{
	switch(typeCase)
	{
	case 1:
		CompleteDataWIGCase1(cell1D, species_set, ND, variabletype, CFD_NT, indexQ, indexV, indexW, indexMIX, indexXs, indexYs);
		break;
	case 2:
		CompleteDataWIGCase2(cell1D, species_set, ND, variabletype, CFD_NT, indexQ, indexV, indexW, indexMIX, indexXs, indexYs);
		break;
	case 3:
		CompleteDataWIGCase3(cell1D, species_set, ND, variabletype, CFD_NT, indexQ, indexV, indexW, indexMIX, indexXs, indexYs);
		break;
	case 4:
		CompleteDataWIGCase4(cell1D, species_set, ND, variabletype, CFD_NT, indexQ, indexV, indexW, indexMIX, indexXs, indexYs);
		break;
	}

	if (is_initialize == true)
	{
		InitializeOtherData(cell1D, CFD_NT, indexR, indexdQ, indexQnew, indexS, indexDivV, indexdpdQ, indexDs, indexMus);
	}
}

void CellData1D::CompleteDataWIG(Data::DataStorageVector<Data::DataStorage>& cell1D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, bool axis, bool viscous, int implicit, int typeCase, bool is_initialize)
{
	int indexQ		= cell1D.dataMap.find(NAME_Q);
	int indexV		= cell1D.dataMap.find(NAME_V);
	int indexW		= cell1D.dataMap.find(NAME_W);
	int indexMIX	= cell1D.dataMap.find(NAME_MIX);
	int indexXs		= cell1D.dataMap.find(NAME_XS);
	int indexYs		= cell1D.dataMap.find(NAME_YS);
	int indexR		= cell1D.dataMap.find(NAME_R);
	int indexdQ		= cell1D.dataMap.find(NAME_dQ);
	int indexQnew	= cell1D.dataMap.find(NAME_Qnew);
	int indexS		= cell1D.dataMap.find(NAME_S);

	int indexDivV	= -1;
	int indexdpdQ	= -1;
	int indexDs		= -1;
	int indexMus	= -1;


	if (axis == true) 	indexDivV	= cell1D.dataMap.find(NAME_DIVV);
	if (implicit != 0)	indexdpdQ	= cell1D.dataMap.find(NAME_dpdQ);
	if (viscous == true)
	{
		indexDs		= cell1D.dataMap.find(NAME_DS);
		indexMus	= cell1D.dataMap.find(NAME_MUS);
	}


	switch(typeCase)
	{
	case 1:
		CompleteDataWIGCase1(cell1D, species_set, ND, variabletype, CFD_NT, indexQ, indexV, indexW, indexMIX, indexXs, indexYs);
		break;
	case 2:
		CompleteDataWIGCase2(cell1D, species_set, ND, variabletype, CFD_NT, indexQ, indexV, indexW, indexMIX, indexXs, indexYs);
		break;
	case 3:
		CompleteDataWIGCase3(cell1D, species_set, ND, variabletype, CFD_NT, indexQ, indexV, indexW, indexMIX, indexXs, indexYs);
		break;
	case 4:
		CompleteDataWIGCase4(cell1D, species_set, ND, variabletype, CFD_NT, indexQ, indexV, indexW, indexMIX, indexXs, indexYs);
		break;
	}

	if (is_initialize == true)
	{
		InitializeOtherData(cell1D, CFD_NT, indexR, indexdQ, indexQnew, indexS, indexDivV, indexdpdQ, indexDs, indexMus);
	}
}



}
}
