/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 15, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_BC_Inviscid_wall.cpp
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


void CFD_BC_Inviscid_wall(CFD_variable_setup_ver2 setup, CELL_CLASS &cell_cl, FACE_CLASS &face, CELL_CLASS &cell_cr, vector<double> &Q_cl, vector<double> &Q_cr)
{
	// 1. Calculate rho*v'
	// Normal rho * vp

	vector<double> rho_V(setup.ND, 0.0);
	vector<double> rho_Vp(setup.ND, 0.0);

	for (int s = 0; s <= setup.ND-1; s++)
	{
		for (int k= 0; k <= setup.ND-1; k++)
		{
			rho_Vp[s]	+= Q_cl[setup.NS+k]*face.n[s][k];
		}
	}

	rho_Vp[0]	= -rho_Vp[0];
	for (int k = 0; k <= setup.ND-1; k++)
	{
		for (int s= 0; s <= setup.ND-1; s++)
		{
			rho_V[k]	+= rho_Vp[s]*face.n[s][k];
		}
	}

	// 2. Apply BC
	for (int s = 0; s <= setup.NS-1; s++)	Q_cr[s]						= Q_cl[s];
	for (int k = 0; k <= setup.ND-1; k++) 	Q_cr[setup.NS+k]		=	 rho_V[k];
	for (int m = 0; m <= setup.NE-1; m++) 	Q_cr[setup.NS+setup.ND+m]	= Q_cl[setup.NS+setup.ND+m];
}



