/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 15, 2015
 *      			Author: Minkwan Kim
 *
 * CalculateJacobianSourceTerm.cpp
 * 			-  
 *  
 */

#include <omp.h>
#include "CFD/include/VariableSet.hpp"
#include "CFD/include/MiscFunctionsCFD.hpp"
#include "CFD/include/FluxInviscid.hpp"
#include "CFD/include/Jacobians.hpp"
#include "Math/include/MathMisc.hpp"

#include "../include/OP2A_Application.hpp"



void ApplicationOP2A::CalculateJacobianSourceTerm()
{
	CalculatedTdQ_dpdQ_at_cells();
	if (problem_setup.is_axisymmetric == true)	JacobianSourceAxis();




	int indexTr;
	int indexTv;
	int indexTe;
	int indexT	= species_set.NS + grid.ND;

	int index_dT = 0;	// Index in cell data2D
	int index_dp = grid.cells[1].data1D.dataMap.find(NAME_dpdQ);
	int VAR = grid.cells[1].data1D(0).numData;

	switch (CFD_variabletype)
	{
	case 1:
		break;

	case 2:
		indexTe = indexT + 1;
		if (problem_setup.is_axisymmetric == true)
		{
#pragma omp parallel for
			for (int c = 1; c <= grid.NCM; c++)
			{
				Data::DataStorage dS("dS_temp", VAR);
				double S = grid.cells[c].geo.S * Math::fabs<double>(grid.cells[c].geo.x[1]);

				CFD::Jacobians::S_he(grid.cells[c].data1D, grid.cells[c].data2D(index_dT), species_set, grid.ND,
						 	 	 	 indexTe, indexQ, indexV, indexW, S, dS);

#pragma ivdep
				for (int i = 0; i <= VAR-1; i++)
				{
					grid.cells[c].data2D(1)(indexTe, i) += dS(i);
				}
			}
		}
		else
		{
#pragma omp parallel for
			for (int c = 1; c <= grid.NCM; c++)
			{
				Data::DataStorage dS("dS_temp", VAR);
				double S = grid.cells[c].geo.S;

				CFD::Jacobians::S_he(grid.cells[c].data1D, grid.cells[c].data2D(index_dT), species_set, grid.ND,
									 indexTe, indexQ, indexV, indexW, S, dS);

#pragma ivdep
				for (int i = 0; i <= VAR-1; i++)
				{
					grid.cells[c].data2D(1)(indexTe, i) += dS(i);
				}
			}
		}
		break;

	case 3:
		indexTv = indexT + 1;

		break;

	case 4:
		indexTv = indexT + 1;
		indexTe = indexT + 2;
		break;

	case 5:
		indexTr = indexT + 1;
		break;

	case 6:
		indexTr = indexT + 1;
		indexTe = indexT + 2;
		break;

	case 7:
		indexTr = indexT + 1;
		indexTv = indexT + 2;
		break;

	case 8:
		indexTr = indexT + 1;
		indexTv = indexT + 2;
		indexTe = indexT + 3;
		break;
	}

}



