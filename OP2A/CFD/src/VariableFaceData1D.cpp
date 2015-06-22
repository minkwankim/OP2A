/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 22, 2015
 *      			Author: Minkwan Kim
 *
 * VariableFaceData1D.cpp
 * 			-  
 *  
 */




#include "CFD/include/Variables.hpp"


namespace OP2A{
namespace CFD{


void FaceData1D::InitializeData(Data::DataStorageVector<Data::DataStorage>& face1D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, bool viscous)
{
	int indexVar;

	indexVar = face1D.dataMap.find(NAME_FINV);
#pragma omp parallel for num_threads(CFD_NT)
	for (int i = 0; i <= face1D(indexVar).numData-1; i++)
		face1D(indexVar).data[i]	= 0.0;


	if (viscous == true)
	{
		indexVar = face1D.dataMap.find(NAME_FVIS);
#pragma omp parallel for num_threads(CFD_NT)
		for (int i = 0; i <= face1D(indexVar).numData-1; i++)
			face1D(indexVar).data[i]	= 0.0;
	}
}







}
}
