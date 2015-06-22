/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 22, 2015
 *      			Author: Minkwan Kim
 *
 * VariableFaceData2D.cpp
 * 			-  
 *  
 */




#include "CFD/include/Variables.hpp"


namespace OP2A{
namespace CFD{


void FaceData2D::InitializeData(Data::DataStorageVector<Data::DataStorage2D>& face2D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, bool viscous, int implicit)
{
	string varname;

	if (implicit != 0)
	{
		int indexVar;
		int numI, numJ;

		// dFinv/dQ +
		varname	= string(NAME_dFINVdQ) + " Plus";
		indexVar = face2D.dataMap.find(varname);
		numI	= face2D.data[indexVar].numData_I;
		numJ	= face2D.data[indexVar].numData_J;

		for (int i = 0; i <= numI-1; i++)
		{
	#pragma omp parallel for num_threads(CFD_NT)
			for (int j = 0; j <= numJ-1; j++)
			{
				face2D.data[indexVar](i,j)	= 0.0;
			}
		}

		// dFinv/dQ -
		varname	= string(NAME_dFINVdQ) + " Minus";
		indexVar = face2D.dataMap.find(varname);
		numI	= face2D.data[indexVar].numData_I;
		numJ	= face2D.data[indexVar].numData_J;

		for (int i = 0; i <= numI-1; i++)
		{
	#pragma omp parallel for num_threads(CFD_NT)
			for (int j = 0; j <= numJ-1; j++)
			{
				face2D.data[indexVar](i,j)	= 0.0;
			}
		}


		if(viscous == true)
		{
			// dFvis/dQ +
			varname = string(NAME_dFVISdQ) + " Plus";
			indexVar = face2D.dataMap.find(varname);
			numI	= face2D.data[indexVar].numData_I;
			numJ	= face2D.data[indexVar].numData_J;

			for (int i = 0; i <= numI-1; i++)
			{
			#pragma omp parallel for num_threads(CFD_NT)
				for (int j = 0; j <= numJ-1; j++)
				{
					face2D.data[indexVar](i,j)	= 0.0;
				}
			}


			// dFvis/dQ -
			varname = string(NAME_dFVISdQ) + " Minus";
			indexVar = face2D.dataMap.find(varname);
			numI	= face2D.data[indexVar].numData_I;
			numJ	= face2D.data[indexVar].numData_J;

			for (int i = 0; i <= numI-1; i++)
			{
			#pragma omp parallel for num_threads(CFD_NT)
				for (int j = 0; j <= numJ-1; j++)
				{
					face2D.data[indexVar](i,j)	= 0.0;
				}
			}
		}
	}
}





}
}
