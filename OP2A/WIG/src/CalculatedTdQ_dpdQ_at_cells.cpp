/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 10, 2015
 *      			Author: Minkwan Kim
 *
 * CalculatedTdQ_dpdQ_at_cells.cpp
 * 			-  
 *  
 */



#include <omp.h>
#include "CFD/include/VariableSet.hpp"
#include "CFD/include/MiscFunctionsCFD.hpp"
#include "CFD/include/FluxInviscid.hpp"

#include "Math/include/MathMisc.hpp"

#include "../include/OP2A_Application.hpp"



void ApplicationOP2A::CalculatedTdQ_dpdQ_at_cells()
{
	int VAR		= grid.cells[1].data1D(indexW).numData;
	int indexE	= species_set.NS + grid.ND;
	int NE		= VAR - indexE;

	int index_dT = 0;	// Index in cell data2D
	int index_dp = grid.cells[1].data1D.dataMap.find(NAME_dpdQ);


#pragma omp parallel for num_threads(CFD_NT)
	for (int c = 1; c <= grid.NCM; c++)
	{
		Data::DataStorageVector<Data::DataStorage>	dT(NE);
		for (int i = 0; i <= NE-1; i++)	dT(i).resize(VAR);

		CFD::Derivatives::dTdQ(grid.cells[c].data1D, species_set, grid.ND, CFD_variabletype, indexQ, indexV, indexW, indexQ, dT);
		CFD::Derivatives::dpdQ(grid.cells[c].data1D, dT, species_set, grid.ND, CFD_variabletype, indexQ, indexV, indexW, indexQ, grid.cells[c].data1D(index_dp));

		for (int e1 = 0; e1 <= NE-1; e1++)
		{
#pragma ivdep
			for (int i = 0; i <= VAR-1; i++)
			{
				grid.cells[c].data2D(index_dT)(e1, i) = dT(e1)(i);
			}
		}
	}
}

