/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 16, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_utilities_residue_norm.cpp
 * 			-  
 *  
 */




#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <omp.h>

#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../../SETUP/include/OP2A_setup_problem.hpp"
#include "../../MATRIX/include/OPPA_Matrix.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../GRID/include/OP2A_grid.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"



void CFD_calc_residue_norms(int NS, int ND, int NE, GRID_CLASS &grid_data, vector< vector<double> > &Residue_c, double &RHS_max, int &RHS_max_i, double &RHS2, int nt)
{
	int 	NCM		= 0;
	int		NCM_p	= grid_data.NCM;

	double 	prhsmax	= -1000.0;
	double 	prhs2	= 0.0;


#pragma omp parallel for reduction(+:prhs2) num_threads(nt)
	for (int c = 0; c <= NCM_p-1; c++)
	{
		double aux_inf = 0.0;
		for (int s = 0; s <= NS-1; s++)	aux_inf += Residue_c[c][s];

		if (prhsmax <= fabs(aux_inf))
		{
			prhsmax 	= fabs(aux_inf);
			RHS_max_i 	= c;
		}

		prhs2 += pow(aux_inf, 2.0);
	}


	// Initialize for non-MPI CASE
	RHS_max		= prhsmax;
	RHS2		= prhs2;
	NCM			= NCM_p;

#ifdef MPI
	MPI_Allreduce(&prhsmax, RHS_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&prhs2, 	   RHS2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&pNCM, 	   &NCM, 1, MPI_INT, 	MPI_SUM, MPI_COMM_WORLD);
#endif

	RHS2 = sqrt((RHS2)/NCM);
}
