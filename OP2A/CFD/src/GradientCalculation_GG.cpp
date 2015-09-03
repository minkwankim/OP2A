/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Sep 1, 2015
 *      			Author: Minkwan Kim
 *
 * GradientCalculation_GG.cpp
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


void GradientCalculation::GradientCalculation_GG(GRID::Grid &grid, int i_VAR1, int i_VAR2, int method, vector<vector <double> >& grad_phi)
{
	int index_count;
	vector<double> phi_f(grid.NFM+1, 0.0);



	switch (method)
	{
	case 0:	// Face Averaging
#pragma omp parallel for
		for (int f = 1; f <= grid.NFM; f++)
		{
			phi_f[f] = (grid.faces[f].geo.cl[0]->data1D.data[i_VAR1].data[i_VAR2]	+  grid.faces[f].geo.cr[0]->data1D.data[i_VAR1].data[i_VAR2]) / 2.0;
		}
		break;


	case 1:	// Inverse Distance Weighted Face Interpolation
#pragma omp parallel for
		for (int f = 1; f <= grid.NFM; f++)
		{
			int N	= grid.faces[f].geo.Neighbor_list.size();
			vector< vector<double> >	di	= Common::vector_2D(N, grid.ND, 0.0);
			vector<double>	di_abs2(N, 0.0);

			for (int i = 0; i <= N-1; i++)
			{
				for (int k = 0; k <= grid.ND-1; k++)
				{
					di[i][k]	= grid.faces[f].geo.Neighbor_list[i]->geo.x[k]	- grid.faces[f].geo.x[k];
					di_abs2[i]	+= di[i][k]*di[i][k];
				}
			}

			double sum_u	= 0.0;
			double sum_b	= 0.0;
			for (int i = 0; i <= N-1; i++)
			{
				sum_u	+= grid.faces[f].geo.Neighbor_list[i]->data1D.data[i_VAR1].data[i_VAR2] / di_abs2[i];
				sum_b	+= 1.0 / di_abs2[i];
			}

			phi_f[f]	= sum_u / sum_b;
		}

		break;


	case 2:	// Weighted Least Squares Face Interpolation
#pragma omp parallel for
		for (int f = 1; f <= grid.NFM; f++)
		{
			int N	= grid.faces[f].geo.Neighbor_list.size();
			vector< vector<double> >	di	= Common::vector_2D(N, grid.ND, 0.0);
			vector<double>	omega_i(N, 0.0);

			for (int i = 0; i <= N-1; i++)
			{
				for (int k = 0; k <= grid.ND-1; k++)
				{
					di[i][k]	= grid.faces[f].geo.Neighbor_list[i]->geo.x[k]	- grid.faces[f].geo.x[k];
					omega_i[i]	+= di[i][k]*di[i][k];
				}

				omega_i[i]	= 1.0 / sqrt(omega_i[i]);
			}

			double sum_omega 		= 0.0;
			double sum_omega_dx 	= 0.0;
			double sum_omega_dy 	= 0.0;
			double sum_omega_dx2 	= 0.0;
			double sum_omega_dy2 	= 0.0;
			double sum_omega_dxdy 	= 0.0;

			for (int i = 0; i <= N-1; i++)
			{
				sum_omega 		+= omega_i[i];
				sum_omega_dx	+= omega_i[i] * di[i][0];
				sum_omega_dy	+= omega_i[i] * di[i][1];
				sum_omega_dx2	+= omega_i[i] * di[i][0]*di[i][0];
				sum_omega_dy2	+= omega_i[i] * di[i][1]*di[i][1];
				sum_omega_dxdy	+= omega_i[i] * di[i][0]*di[i][1];

			}


			Math::MATRIX A(grid.ND+1, grid.ND+1, false);

			if (grid.ND == 2)
			{
				A(0,0)	= sum_omega;
				A(0,1)	= sum_omega_dx;
				A(0,2)	= sum_omega_dy;

				A(1,0)	= sum_omega_dx;
				A(1,1)	= sum_omega_dx2;
				A(1,2)	= sum_omega_dxdy;

				A(2,0)	= sum_omega_dy;
				A(2,1)	= sum_omega_dxdy;
				A(2,2)	= sum_omega_dy2;
			}
			else
			{
				double sum_omega_dz 	= 0.0;
				double sum_omega_dz2 	= 0.0;
				double sum_omega_dxdz 	= 0.0;
				double sum_omega_dydz 	= 0.0;

				for (int i = 0; i <= N-1; i++)
				{
					sum_omega_dz	+= omega_i[i] * di[i][2];
					sum_omega_dz2	+= omega_i[i] * di[i][2]*di[i][2];
					sum_omega_dxdz	+= omega_i[i] * di[i][0]*di[i][2];
					sum_omega_dydz	+= omega_i[i] * di[i][1]*di[i][2];
				}

				A(0,0)	= sum_omega;
				A(0,1)	= sum_omega_dx;
				A(0,2)	= sum_omega_dy;
				A(0,3)	= sum_omega_dz;

				A(1,0)	= sum_omega_dx;
				A(1,1)	= sum_omega_dx2;
				A(1,2)	= sum_omega_dxdy;
				A(1,3)  = sum_omega_dxdz;

				A(2,0)	= sum_omega_dy;
				A(2,1)	= sum_omega_dxdy;
				A(2,2)	= sum_omega_dy2;
				A(2,3)  = sum_omega_dydz;

				A(3,0)	= sum_omega_dz;
				A(3,1)	= sum_omega_dxdz;
				A(3,2)	= sum_omega_dydz;
				A(3,3)	= sum_omega_dz2;
			}


			Math::MATRIX Ainv = Math::MATRIX_Inv(A);
			Math::MATRIX B(grid.ND+1, N, false);

			for (int i = 0; i <= N-1; i++)
			{
				B(0, i) = omega_i[i];
				for (int k = 0; k <= grid.ND-1; k++) B(1+k, i) = omega_i[i]*di[i][k];
			}

			vector<double> C(N, 0.0);

#pragma ivdep
			for (int i = 0; i <= N-1; i++)
			{
				for (int k = 0; k <= grid.ND; k++)
				{
					C[i]	+= Ainv(0,k)*B(k,i);
				}
			}

			for (int i = 0; i <= N-1; i++)
			{
				phi_f[f] += C[i] * grid.faces[f].geo.Neighbor_list[i]->data1D.data[i_VAR1].data[i_VAR2];
			}
		}

		break;
	}


	for (int c = 1; c <= grid.NCM; c++)
	{
		for (int k = 0; k <= grid.ND-1; k++)
		{
			grad_phi[c][k]	= 0.0;
			for (int f = 0; f <= grid.cells[c].geo.NF-1; f++)
			{
				if (grid.cells[c].geo.ID == grid.cells[c].geo.face_list[f]->geo.cl[0]->geo.ID)
				{
					grad_phi[c][k]	+= phi_f[grid.cells[c].geo.face_list[f]->geo.ID] * grid.cells[c].geo.face_list[f]->geo.n[0][k] * grid.cells[c].geo.face_list[f]->geo.S;
				}
				else
				{
					grad_phi[c][k]	-= phi_f[grid.cells[c].geo.face_list[f]->geo.ID] * grid.cells[c].geo.face_list[f]->geo.n[0][k] * grid.cells[c].geo.face_list[f]->geo.S;
				}
			}

			grad_phi[c][k]	/= grid.cells[c].geo.S;
		}
	}




}

}
}
