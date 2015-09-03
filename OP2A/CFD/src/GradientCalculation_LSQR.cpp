/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Sep 1, 2015
 *      			Author: Minkwan Kim
 *
 * GradientCalculation_LSQR.cpp
 * 			-  
 *  
 */



#include <omp.h>
#include <time.h>

#include "CFD/include/GradientCalculation.hpp"

#include "Math/include/OP2A_Matrix.hpp"
#include "Math/include/MathMisc.hpp"
#include "Common/include/Time_Info.hpp"
#include "Common/include/MultiDimension.hpp"

namespace OP2A{
namespace CFD{


void GradientCalculation::GradientCalculation_LSQR(GRID::Grid &grid, int i_VAR1, int i_VAR2, int method, vector<vector <double> >& grad_phi)
{
#pragma omp parallel for
	for (int c = 1; c <= grid.NCM; c++)
	{
		int N = grid.cells[c].geo.neighbor_list.size();
		vector<double>	B(grid.ND, 0.0);

		for (int i = 0; i <= N-1; i++)
		{
			double delta_phi	= grid.cells[c].geo.neighbor_list[i]->data1D.data[i_VAR1].data[i_VAR2] - grid.cells[c].data1D.data[i_VAR1].data[i_VAR2];

			for (int k = 0; k <= grid.ND-1; k++)
			{
				B[k]	+= grid.cells[c].geo.omega_i2_dx[i][k] * delta_phi;
			}
		}


		for (int k = 0; k <= grid.ND-1; k++)
		{
			grad_phi[c][k]	= 0.0;
			for (int k1 = 0; k1 <= grid.ND-1; k1++)
			{
				grad_phi[c][k]	+= grid.cells[c].geo.LSQRS_matrix[k][k1]*B[k1];
			}
		}
	}
}




}
}
