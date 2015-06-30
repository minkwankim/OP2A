/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 29, 2015
 *      			Author: Minkwan Kim
 *
 * TimeIntegrateExplicit.cpp
 * 			-  
 *  
 */



#include <omp.h>
#include "CFD/include/VariableSet.hpp"
#include "CFD/include/MiscFunctionsCFD.hpp"
#include "CFD/include/FluxInviscid.hpp"

#include "Math/include/MathMisc.hpp"

#include "../include/OP2A_Application.hpp"


void ApplicationOP2A::TimeIntegrateExplicit()
{
	if (problem_setup.is_axisymmetric == true)
	{
#pragma omp parallel for num_threads(CFD_NT)
		for (int c = 1; c <= grid.NCM; c++)
		{
			double Vol 		= grid.cells[c].geo.S * Math::fabs<double>(grid.cells[c].geo.x[1]);
			double dt_Vol	= dt / Vol;

			for (int i = 0; i <= grid.cells[c].data1D(indexResidue).numData-1; i++)
			{
				grid.cells[c].data1D(indexdQ)(i) = -dt_Vol * grid.cells[c].data1D(indexResidue)(i);
			}
		}
	}
	else
	{
#pragma omp parallel for num_threads(CFD_NT)
		for (int c = 1; c <= grid.NCM; c++)
		{
			double Vol 		= grid.cells[c].geo.S;
			double dt_Vol	= dt / Vol;

			for (int i = 0; i <= grid.cells[c].data1D(indexResidue).numData-1; i++)
			{
				grid.cells[c].data1D(indexdQ)(i) = -dt_Vol * grid.cells[c].data1D(indexResidue)(i);
			}
		}
	}
}
