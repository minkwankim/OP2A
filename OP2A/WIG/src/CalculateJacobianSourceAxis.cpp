/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 10, 2015
 *      			Author: Minkwan Kim
 *
 * CalculateJacobianSourceAxis.cpp
 * 			-  
 *  
 */


#include <omp.h>
#include "CFD/include/VariableSet.hpp"
#include "CFD/include/MiscFunctionsCFD.hpp"
#include "CFD/include/Jacobians.hpp"

#include "Math/include/MathMisc.hpp"
#include "../include/OP2A_Application.hpp"


void ApplicationOP2A::JacobianSourceAxis()
{
	int indexdSdQ = 1;
	int indexdpdQ = grid.cells[1].data1D.dataMap.find(NAME_dpdQ);

	int NS	= species_set.NS;
	int NE	= grid.cells[1].data1D(0).numData - NS - grid.ND;

#pragma omp parallel for num_threads(CFD_NT)
	for (int c = 1; c <= grid.NCM; c++)
	{
		CFD::Jacobians::Axisymmertic(grid.cells[c].data1D(indexdpdQ), NS, grid.ND, NE, grid.cells[c].geo.S, grid.cells[c].data2D(indexdSdQ));
	}
}

