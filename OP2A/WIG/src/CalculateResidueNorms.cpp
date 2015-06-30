/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 29, 2015
 *      			Author: Minkwan Kim
 *
 * CalculateResidueNorms.cpp
 * 			-  
 *  
 */


#include <omp.h>
#include "CFD/include/VariableSet.hpp"
#include "CFD/include/MiscFunctionsCFD.hpp"
#include "../include/OP2A_Application.hpp"
#include "Math/include/MathMisc.hpp"



void ApplicationOP2A::CalcualtedResidueNorms()
{
	double 	prhsmax	= -1000.0;
	double 	prhs2	= 0.0;

#pragma omp parallel for reduction(+:prhs2)
	for (int c = 1; c <= grid.NCM; c++)
	{
		double aux_inf = 0.0;
#pragma omp parallel for reduction(+:aux_inf)
		for (int s = 0; s <= species_set.NS-1; s++)
		{
			aux_inf += grid.cells[c].data1D(indexResidue)(s);
		}

		if (prhsmax <= fabs(aux_inf))
		{
			prhsmax = Math::fabs<double>(aux_inf);
			RHS_n 	= c;
		}

		prhs2 += pow(aux_inf, 2.0);
	}


	// Initialize for non-MPI CASE
	RHS_max		= prhsmax;
	RHS_2		= prhs2;

#ifdef MPI
	MPI_Allreduce(&prhsmax, RHS_max,	1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&prhs2,	RHS_2, 		1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&grid.NCM, &grid.NCM, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

	RHS_2 = sqrt((RHS_2)/grid.NCM);


	if (RHS_max <= problem_setup.convergence_criterion && problem_setup.n_current >= 100)
	{
		termination = true;
	}
	else
	{
		termination = false;
	}
}
