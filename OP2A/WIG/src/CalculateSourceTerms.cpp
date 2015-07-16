/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 15, 2015
 *      			Author: Minkwan Kim
 *
 * CalculateSourceTerms.cpp
 * 			-  
 *  
 */



#include <omp.h>
#include <limits>
#include "Common/include/Exception_NaNValue.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"
#include "CFD/include/Source_NONEQ.hpp"

#include "Math/include/MathMisc.hpp"
#include "../include/OP2A_Application.hpp"

void ApplicationOP2A::CalculateSourceTerm()
{
	// Calcualete course term for Thermal NONEQ
	CalculateSourceNONEQ(true);





	// Update Residue
	if (problem_setup.is_axisymmetric == true)
	{
#pragma omp parallel for num_threads(CFD_NT)
		for (int c = 1; c <= grid.NCM; c++)
		{
			double S = 0.0;
			S = grid.cells[c].geo.S * Math::fabs<double>(grid.cells[c].geo.x[1]);

#pragma ivdep
			for (int s = 0; s <= grid.cells[c].data1D(indexResidue).numData-1; s++)
			{
				grid.cells[c].data1D(indexResidue)(s) += -S*grid.cells[c].data1D(indexS)(s);
			}

		}
	}
	else
	{
#pragma omp parallel for num_threads(CFD_NT)
		for (int c = 1; c <= grid.NCM; c++)
		{
			double S = 0.0;
			S = grid.cells[c].geo.S;

#pragma ivdep
			for (int s = 0; s <= grid.cells[c].data1D(indexResidue).numData-1; s++)
			{
				grid.cells[c].data1D(indexResidue)(s) += -S*grid.cells[c].data1D(indexS)(s);
			}
		}
	}
}
