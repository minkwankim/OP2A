/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 15, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_BC_Inviscid_inlet.cpp
 * 			-  
 *  
 */


#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include "../include/OP2A_CFD_Variable_change.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../../GRID/include/OP2A_grid.hpp"


void CFD_BC_Inviscid_inlet(CFD_variable_setup_ver2 setup, CELL_CLASS &cell_cl, CELL_CLASS &cell_cr, vector<double> &Q_cl, vector<double> &Q_cr, vector<double> &Q_inlet)
{
	// 1. Apply BC
#pragma omp parallel for
	for (int i = 0; i <= setup.VAR-1; i++)	Q_cr[i]			= Q_inlet[i];
}
