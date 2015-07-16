/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 9, 2015
 *      			Author: Minkwan Kim
 *
 * UpdateFluxJacobian.cpp
 * 			-  
 *  
 */



#include <omp.h>
#include "CFD/include/VariableSet.hpp"
#include "CFD/include/MiscFunctionsCFD.hpp"
#include "CFD/include/FluxInviscid.hpp"

#include "Math/include/MathMisc.hpp"

#include "../include/OP2A_Application.hpp"


void ApplicationOP2A::UpdateFluxJacobian()
{
	int VAR = grid.cells[1].data1D(0).numData;

	if (problem_setup.is_axisymmetric == true)
	{
#pragma omp parallel for num_threads(CFD_NT)
		for (int f = 1; f <= grid.NFM; f++)
		{
			double S = grid.faces[f].geo.S * Math::fabs<double>(grid.faces[f].geo.x[1]);

#pragma ivdep
			for (int i = 0; i <= VAR-1; i++)
			{
				for (int j = 0; j <= VAR-1; j++)
				{
					grid.faces[f].data2D(0)(i,j) *= S;
					grid.faces[f].data2D(1)(i,j) *= S;
				}
			}
		}
	}
	else
	{
#pragma omp parallel for num_threads(CFD_NT)
		for (int f = 1; f <= grid.NFM; f++)
		{
			double S = grid.faces[f].geo.S;

#pragma ivdep
			for (int i = 0; i <= VAR-1; i++)
			{
				for (int j = 0; j <= VAR-1; j++)
				{
					grid.faces[f].data2D(0)(i,j) *= S;
					grid.faces[f].data2D(1)(i,j) *= S;
				}
			}
		}
	}
}
