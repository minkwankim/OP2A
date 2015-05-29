/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 15, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_utilities_time.cpp
 * 			-  
 *  
 */



#include <mkl.h>
#include <vector>
#include <limits>
#include "omp.h"
#include "../../MATRIX/include/math_misc.hpp"

using namespace std;

double calculate_dt(int NCM, vector<double> &dX, vector<double> &U, double dt_MAX, double CFL)
{
	double dt;
	double min_dt		= 1.0e99;
	double min_dt_mpi	= min_dt;



	for (int i = 0; i <= NCM-1; i++)
	{
		dt		= CFL * dX[i] / U[i];
		min_dt	= fmin(dt, min_dt);
	}

	if (min_dt > dt_MAX)	min_dt	= dt_MAX;
	if (min_dt < 1.0e-14)	min_dt	= 1.0e-14;


#ifdef MPI
	min_dt_mpi	= min_dt;
	MPI_Allreduce(&min_dt_mpi, &min_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

	return (min_dt);
}



double calculate_dt_multi(int NCM, vector<double> &dX, vector< vector<double> > &U, double max_dt, double CFL, int NUM_FLUID)
{
	double dt	= 0.0;
	vector<double>				dt_temp(NUM_FLUID, 0.0);

#pragma omp parallel for num_threads(NUM_FLUID)
	for (int f = 0; f <= NUM_FLUID-1; f++)
		dt_temp[f] = calculate_dt(NCM, dX, U[f], max_dt, CFL);

	dt = dt_temp[0];
	for (int f = 1; f <= NUM_FLUID-1; f++)	dt = fmin(dt, dt_temp[f]);

	return (dt);
}
