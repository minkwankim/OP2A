/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 16, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Time_integral_Implicit_point.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_CFD_Variable_change.hpp"
#include "../include/OP2A_CFD_reconstruct.hpp"
#include "../include/OP2A_CFD_Flux.hpp"


#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../../SETUP/include/OP2A_setup_problem.hpp"
#include "../../MATRIX/include/OPPA_Matrix.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../GRID/include/OP2A_grid.hpp"
#include "../../MATRIX/include/CMatrix.hpp"

void CFD_time_integration_Implicit_point(GRID_CLASS &grid, SOL_CFD &Solution, double dt, bool is_axis, int nt)
{
	int NVAR	= Solution.setup.VAR;
	int NCM_NGM = grid.NCM + grid.NGM;

	double flag_axis;
	if (is_axis	== true)	flag_axis	= 1.0;
	else					flag_axis	= 0.0;


	// Initialize Matrix
	vector<CMatrix_ver2>	M_cl(grid.NCM);
	vector<CMatrix_ver2>	B(grid.NCM);
	vector<CMatrix_ver2>	R(grid.NCM);
	vector<CMatrix_ver2>	X(grid.NCM);

	for (int i = 0; i <= grid.NCM-1; i++)
	{
		M_cl[i].zeros(NVAR, NVAR);
		B[i].zeros(NVAR, NVAR);

		R[i].zeros(NVAR, 1);
		X[i].zeros(NVAR, 1);
	}




	// STEP 1: Set-up matrix M_CL and R_CL
	// [!!!! NOTE !!!!] Index range of CMatrix is same as Matlab. It starts with 1 (Not ZERO)
#pragma omp parallel for num_threads(nt)
	for (int i = 0; i <= grid.NCM-1; i++)
	{
		double Volume	= grid.cells.data_ptr[i]->S * (1.0 + flag_axis*(grid.cells.data_ptr[i]->x[1]- 1.0));
		double Vol_dt	= Volume / dt;

		for (int r = 1; r <= NVAR; r++)
			for (int l = 1; l <= NVAR; l++) M_cl[i](r, l)	= -Solution.Jacobian_source[i][r-1][l-1];

		for (int r = 1; r <= NVAR; r++)	M_cl[i](r, r)	+= Vol_dt;
		for (int r = 1; r <= NVAR; r++)	R[i](r,1)		= -Solution.Rn[i][r-1];
	}


	// SETP 2: ADD Viscous/inviscid Jacobian
#pragma omp parallel for num_threads(nt)
	for (int c = 1; c <= grid.NCM; c++)
	{
		int c_ptr	= grid.cells.whereis[c];
		int cl, cr;
		int f;

		for (int k = 0; k <= grid.cells.data_ptr[c_ptr]->NF-1; k++)
		{
			f 	= grid.faces.whereis[grid.cells.data_ptr[c_ptr]->face[k]];
			cl	= grid.faces.data_ptr[f]->cl[0];
			cr	= grid.faces.data_ptr[f]->cr[0];

			if (cl == c)
			{
				for (int j = 0; j <= NVAR-1; j++)
					for (int k = 0; k <= NVAR-1; k++)
						M_cl[c_ptr](j+1, k+1) += Solution.Jacobian_inviscid_plus[f][j][k];
			}

			if (cr == c)
			{
				for (int j = 0; j <= NVAR-1; j++)
					for (int k = 0; k <= NVAR-1; k++)
						M_cl[c_ptr](j+1, k+1) += -Solution.Jacobian_inviscid_minus[f][j][k];;
			}
		}
	}








	// STEP 6 SOLVE BLOCK TRI-DIAGONAL MATRIX
	// 6.1. Initialize
	B 	= M_cl;
	vector<CMatrix_ver2>	B_inv(grid.NCM);

#pragma omp parallel for num_threads(nt)
	for (int i = 0; i <= grid.NCM-1; i++)	B_inv[i]	= CMatrix_INV(B[i]);


	// 6.2 Solve
	// 		1. Initial solution
#pragma omp parallel for num_threads(nt)
	for (int i = 0; i <= grid.NCM-1; i++)	X[i] = B_inv[i] * R[i];

#pragma omp parallel for num_threads(nt)
	for (int i = 0; i <= grid.NCM-1; i++)
		for (int j = 0; j <= NVAR-1; j++)
			Solution.dQ[i][j]	= X[i](j+1, 1);


	for (int p = 1; p <= 4; p++)
	{
#pragma omp parallel for num_threads(nt)
		for (int i = 0; i <= grid.NCM-1; i++)
			for (int r = 0; r <= NVAR-1; r++)	R[i](r+1,1) = -Solution.Rn[i][r];


		// STep 2 Communication (FOR MPI)


		// Step 3. Update R(n+1, p)
#pragma omp parallel for num_threads(nt)
		for (int c = 1; c <= grid.NCM; c++)
		{
			int c_ptr	= grid.cells.whereis[c];
			int cl, cr;
			int f;

			for (int k = 0; k <= grid.cells.data_ptr[c_ptr]->NF-1; k++)
			{
				f 	= grid.faces.whereis[grid.cells.data_ptr[c_ptr]->face[k]];
				cl	= grid.faces.data_ptr[f]->cl[0];
				cr	= grid.faces.data_ptr[f]->cr[0];
				int c_ptr2;

				if (cl == c)
				{
					if (cr > 0)
					{
						c_ptr2	= grid.cells.whereis[cr];
						for (int j = 0; j <= NVAR-1; j++)
						{
							double aux = 0.0;
							for (int l = 0; l <= NVAR-1; l++)
								aux +=	Solution.Jacobian_inviscid_minus[f][j][l] * Solution.dQ[c_ptr2][l];

							R[c_ptr](j+1, 1)	-= aux;
						}
					}
				}

				if (cr == c)
				{
					if (cl > 0)
					{
						c_ptr2	= grid.cells.whereis[cl];
						for (int j = 0; j <= NVAR-1; j++)
						{
							double aux = 0.0;
							for (int l = 0; l <= NVAR-1; l++)
								aux +=	Solution.Jacobian_inviscid_plus[f][j][l] * Solution.dQ[c_ptr2][l];

							R[c_ptr](j+1, 1)	+= aux;
						}
					}
				}
			}
		}

#pragma omp parallel for num_threads(nt)
		for (int i = 0; i <= grid.NCM-1; i++)	X[i] = B_inv[i] * R[i];

#pragma omp parallel for num_threads(nt)
		for (int i = 0; i <= grid.NCM-1; i++)
			for (int r = 0; r <= NVAR-1; r++)
				Solution.dQ[i][r]	= X[i](r+1, 1);
	}

}
