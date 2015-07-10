/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 9, 2015
 *      			Author: Minkwan Kim
 *
 * TimeIntegrateImplicitPoint.cpp
 * 			-  
 *  
 */


#include <omp.h>
#include "CFD/include/VariableSet.hpp"
#include "CFD/include/MiscFunctionsCFD.hpp"
#include "CFD/include/FluxInviscid.hpp"

#include "Math/include/MathMisc.hpp"
#include "Math/include/OP2A_Matrix.hpp"
#include "../include/OP2A_Application.hpp"


void ApplicationOP2A::TimeIntegrateImplicitPoint()
{
	int VAR	= grid.cells[1].data1D(0).numData;
	int index_J_inv_plus  = 0;
	int index_J_inv_minus = 1;


	// [PRE] Initialize Matrix
	Math::MATRIX	matrix_temp1(VAR, VAR, false);	matrix_temp1.zeros();
	Math::MATRIX	matrix_temp2(VAR, 1, false);	matrix_temp2.zeros();

	vector<Math::MATRIX> M_cl(grid.NCM+1, matrix_temp1);
	vector<Math::MATRIX> B(grid.NCM+1, matrix_temp1);

	vector<Math::MATRIX> R(grid.NCM+1, matrix_temp2);
	vector<Math::MATRIX> X(grid.NCM+1, matrix_temp2);


	// STEP 1: Set-up matrix M_CL and R_CL
	if (problem_setup.is_axisymmetric == true)
	{
#pragma omp parallel for num_threads(CFD_NT)
		for (int c = 1; c <= grid.NCM; c++)
		{
			double Vol 		= grid.cells[c].geo.S * Math::fabs<double>(grid.cells[c].geo.x[1]);
			double Vol_dt	= Vol / dt;

			for (int r = 0; r <= VAR-1; r++)
			{
				for (int l = 0; l <= VAR-1; l++)
				{
					M_cl[c](r, l)	= -grid.cells[c].data2D(1)(r,l);
				}
			}

			for (int r = 0; r <= VAR-1; r++)
			{
				M_cl[c](r, r)	+= Vol_dt;
			}

			for (int r = 0; r <= VAR-1; r++)
			{
				R[c](r,0) = -grid.cells[c].data1D(indexResidue)(r);
			}
		}
	}
	else
	{
#pragma omp parallel for num_threads(CFD_NT)
		for (int c = 1; c <= grid.NCM; c++)
		{
			double Vol 		= grid.cells[c].geo.S;
			double Vol_dt	= Vol / dt;

			for (int r = 0; r <= VAR-1; r++)
			{
				for (int l = 0; l <= VAR-1; l++)
				{
					M_cl[c](r, l)	= -grid.cells[c].data2D(1)(r,l);
				}
			}

			for (int r = 0; r <= VAR-1; r++)
			{
				M_cl[c](r, r)	+= Vol_dt;
			}

			for (int r = 0; r <= VAR-1; r++)
			{
				R[c](r,0) = -grid.cells[c].data1D(indexResidue)(r);
			}
		}
	}



	// SETP 2: ADD Viscous/inviscid Jacobian
#pragma omp parallel for num_threads(CFD_NT)
	for (int c = 1; c <= grid.NCM; c++)
	{
		for (int f = 0; f <= grid.cells[c].geo.NF-1; f++)
		{
			if (grid.cells[c].geo.face_list[f]->geo.cl[0]->geo.ID == grid.cells[c].geo.ID)
			{
				for (int j = 0; j <= VAR-1; j++)
				{
					for (int k = 0; k <= VAR-1; k++)
					{
						M_cl[c](j, k) += grid.cells[c].geo.face_list[f]->data2D(index_J_inv_plus)(j,k);
					}
				}
			}


			if (grid.cells[c].geo.face_list[f]->geo.cr[0]->geo.ID == grid.cells[c].geo.ID)
			{
				for (int j = 0; j <= VAR-1; j++)
				{
					for (int k = 0; k <= VAR-1; k++)
					{
						M_cl[c](j, k) += -grid.cells[c].geo.face_list[f]->data2D(index_J_inv_minus)(j,k);
					}
				}
			}
		}

	}


	/*
	 * @todo: Need to add Viscous
	 */










	// STEP 6 SOLVE BLOCK TRI-DIAGONAL MATRIX
	// 6.1. Initialize
	B 	= M_cl;
	vector<Math::MATRIX> B_inv(grid.NCM+1, matrix_temp1);

#pragma omp parallel for num_threads(CFD_NT)
	for (int i = 1; i <= grid.NCM; i++)
	{
		B_inv[i]	= MATRIX_Inv(B[i]);
	}

	// 6.2 Solve
	// 		1. Initial solution
#pragma omp parallel for num_threads(CFD_NT)
	for (int i = 1; i <= grid.NCM; i++)
	{
		X[i] = B_inv[i] * R[i];
		for (int r = 0; r <= VAR-1; r++)
		{
			if(X[i](r,0) != X[i](r,0))
			{
				std::ostringstream oss;
				oss << "It has NaN value: [Cell ID]: " << i << "  [VAR]: " << r;
				throw Common::ExceptionNaNValue (FromHere(), oss.str());
			}
		}

	}

#pragma omp parallel for num_threads(CFD_NT)
	for (int c = 1; c <= grid.NCM; c++)
	{
		for (int j = 0; j <= VAR-1; j++)
		{
			grid.cells[c].data1D(indexdQ)(j) = X[c](j, 0);
		}
	}




	for (int p = 1; p <= 4; p++)
	{
#pragma omp parallel for num_threads(CFD_NT)
		for (int c = 1; c <= grid.NCM; c++)
		{
			for (int r = 0; r <= VAR-1; r++)
			{
				R[c](r,0) = -grid.cells[c].data1D(indexResidue)(r);
			}
		}


		// STep 2 Communication (FOR MPI)
		// @todo: Need to update for MPI



		// Step 3. Update R(n, p)
#pragma omp parallel for num_threads(CFD_NT)
		for (int c = 1; c <= grid.NCM; c++)
		{
			for (int k = 0; k <= grid.cells[c].geo.NF-1; k++)
			{
				if (grid.cells[c].geo.face_list[k]->geo.cl[0]->geo.ID == grid.cells[c].geo.ID)
				{
					if (grid.cells[c].geo.face_list[k]->geo.cr[0]->geo.ID > 0)
					{
						for (int j = 0; j <= VAR-1; j++)
						{
							double aux = 0.0;
							for (int l = 0; l <= VAR-1; l++)
							{
								aux +=	grid.cells[c].geo.face_list[k]->data2D(index_J_inv_minus)(j,l) * grid.cells[c].geo.face_list[k]->geo.cr[0]->data1D(indexdQ)(l);
							}

							R[c](j, 0)	-= aux;
						}
					}
				}


				if (grid.cells[c].geo.face_list[k]->geo.cr[0]->geo.ID == grid.cells[c].geo.ID)
				{
					if (grid.cells[c].geo.face_list[k]->geo.cl[0]->geo.ID > 0)
					{
						for (int j = 0; j <= VAR-1; j++)
						{
							double aux = 0.0;
							for (int l = 0; l <= VAR-1; l++)
							{
								aux +=	grid.cells[c].geo.face_list[k]->data2D(index_J_inv_plus)(j,l) * grid.cells[c].geo.face_list[k]->geo.cl[0]->data1D(indexdQ)(l);
							}

							R[c](j, 0)	+= aux;
						}
					}
				}
			}
		}

#pragma omp parallel for num_threads(CFD_NT)
		for (int i = 1; i <= grid.NCM; i++)
		{
			X[i] = B_inv[i] * R[i];
		}

#pragma omp parallel for num_threads(CFD_NT)
		for (int i = 1; i <= grid.NCM; i++)
		{
			for (int r = 0; r <= VAR-1; r++)
			{
				if(X[i](r,0) != X[i](r,0))
				{
					std::ostringstream oss;
					oss << "It has NaN value: [Cell ID]: " << i << "  [VAR]: " << r;
					throw Common::ExceptionNaNValue (FromHere(), oss.str());
				}

				grid.cells[i].data1D(indexdQ)(r) = X[i](r,0);
			}
		}
	}
}
