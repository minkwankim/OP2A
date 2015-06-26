/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 24, 2015
 *      			Author: Minkwan Kim
 *
 * Calculatedt.cpp
 * 			-  
 *  
 */



#include <omp.h>
#include "CFD/include/VariableSet.hpp"
#include "CFD/include/MiscFunctionsCFD.hpp"
#include "../include/OP2A_Application.hpp"
#include "Math/include/MathMisc.hpp"


void ApplicationOP2A::Calcualtedt()
{
	int index_u	= grid.cells[1].data1D(indexV).dataMap.find(CFD::CFD_VariableSet::var_stringV(0, CFD::FluxCategory::Momentum));
	double dt_min	= 1.0e-14;

	vector<double>	dt_temp(grid.NCM+1, dt_min);

#pragma omp parallel for
	for (int c = 1; c <= grid.NCM; c++)
	{
		double U_abs = 0.0;
		for (int k = 0; k <= grid.ND-1; k++)
		{
			U_abs	+= pow(grid.cells[c].data1D(indexV)(index_u+k), 2.0);
		}
		U_abs = sqrt(U_abs);

		dt_temp[c]	= CFD::Calculatedt(grid.cells[c].geo.characteristic_length, U_abs, problem_setup.max_dt, CFLNumber);
	}

	dt	= Math::fminOMP(dt_temp, 1, grid.NCM);
}
