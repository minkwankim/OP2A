/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 22, 2015
 *      			Author: Minkwan Kim
 *
 * VariableCellData2D.cpp
 * 			-  
 *  
 */




#include "CFD/include/Variables.hpp"


namespace OP2A{
namespace CFD{


void CellData2D::InitializeData(Data::DataStorageVector<Data::DataStorage2D>& cell2D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, bool viscous, int implicit)
{
	int indexVar;
	int numI;
	int numJ;

	if (viscous == true)
	{
		indexVar = cell2D.dataMap.find(NAME_KAPPAS);
		numI	= cell2D.data[indexVar].numData_I;
		numJ	= cell2D.data[indexVar].numData_J;

		for (int i = 0; i <= numI-1; i++)
		{
#pragma omp parallel for num_threads(CFD_NT)
			for (int j = 0; j <= numJ-1; j++)
			{
				cell2D.data[indexVar](i,j)	= 0.0;
			}
		}
	}

	if (implicit != 0)
	{
		indexVar = cell2D.dataMap.find(NAME_dTdQ);
		numI	= cell2D.data[indexVar].numData_I;
		numJ	= cell2D.data[indexVar].numData_J;

		for (int i = 0; i <= numI-1; i++)
		{
#pragma omp parallel for num_threads(CFD_NT)
			for (int j = 0; j <= numJ-1; j++)
			{
				cell2D.data[indexVar](i,j)	= 0.0;
			}
		}

		indexVar = cell2D.dataMap.find(NAME_dSdQ);
		numI	= cell2D.data[indexVar].numData_I;
		numJ	= cell2D.data[indexVar].numData_J;

		for (int i = 0; i <= numI-1; i++)
		{
#pragma omp parallel for num_threads(CFD_NT)
			for (int j = 0; j <= numJ-1; j++)
			{
				cell2D.data[indexVar](i,j)	= 0.0;
			}
		}
	}

}


void CellData2D::InitializeData(Data::DataStorageVector<Data::DataStorage2D>& cell2D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, int indexdTdQ, int indexdSdQ, int indexKappas)
{
	int numI;
	int numJ;

	if (indexdTdQ >= 0)
	{
		numI	= cell2D.data[indexdTdQ].numData_I;
		numJ	= cell2D.data[indexdTdQ].numData_J;

		for (int i = 0; i <= numI-1; i++)
		{
#pragma omp parallel for num_threads(CFD_NT)
			for (int j = 0; j <= numJ-1; j++)
			{
				cell2D.data[indexdTdQ](i,j)	= 0.0;
			}
		}
	}


	if (indexdSdQ >= 0)
	{
		numI	= cell2D.data[indexdSdQ].numData_I;
		numJ	= cell2D.data[indexdSdQ].numData_J;

		for (int i = 0; i <= numI-1; i++)
		{
#pragma omp parallel for num_threads(CFD_NT)
			for (int j = 0; j <= numJ-1; j++)
			{
				cell2D.data[indexdSdQ](i,j)	= 0.0;
			}
		}
	}



	if (indexKappas >= 0)
	{
		numI	= cell2D.data[indexKappas].numData_I;
		numJ	= cell2D.data[indexKappas].numData_J;

		for (int i = 0; i <= numI-1; i++)
		{
#pragma omp parallel for num_threads(CFD_NT)
			for (int j = 0; j <= numJ-1; j++)
			{
				cell2D.data[indexKappas](i,j)	= 0.0;
			}
		}
	}

}



}
}
