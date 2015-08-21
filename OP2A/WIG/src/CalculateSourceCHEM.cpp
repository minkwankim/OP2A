/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 17, 2015
 *      			Author: Minkwan Kim
 *
 * CalculateSourceCHEM.cpp
 * 			-  
 *  
 */



#include <omp.h>
#include <limits>
#include "Common/include/Exception_NaNValue.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"
#include "CFD/include/Source_CHEM.hpp"

#include "Math/include/MathMisc.hpp"
#include "../include/OP2A_Application.hpp"


void ApplicationOP2A::CalculateSourceCHEM(bool needToInitialize)
{
	if (needToInitialize == true)
	{
#pragma omp parallel
		for (int c = 1; c <= grid.NCM; c++)
		{
	#pragma ivdep
			for (int i = 0; i <= grid.cells[c].data1D(indexS).numData-1; i++)
			{
				grid.cells[c].data1D(indexS)(i)	= 0.0;
			}
		}
	}


	// Calculate CHEM Source term
#pragma omp parallel for num_threads(CFD_NT)
	for (int c = 1; c <= grid.NCM; c++)
	{
		CFD::Source_CHEM::S_chem(grid.cells[c].data1D, species_set, grid.ND, CFD_variabletype, indexQ, indexV, indexW, grid.cells[c].data1D(indexS));
	}
}
