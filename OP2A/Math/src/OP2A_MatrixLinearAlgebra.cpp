/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 26, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_MatrixLinearAlgebra.cpp
 * 			-  
 *  
 */

#include "mkl.h"
#include "Math/include/OP2A_Matrix.hpp"
#include "Common/include/Exception_MemoryAllocator.hpp"
#include "Common/include/Exception_DimensionMatch.hpp"


namespace OP2A{
namespace Math{

MATRIX& MATRIX::MINOR(const int r, const int c)
{
	removeRow(r);
	removeColumn(c);

	return *this;
}

MATRIX MATRIX_Minor(MATRIX &A, const int r, const int c)
{
	MATRIX C = A;
	C.MINOR(r, c);
	return (C);
}



double MATRIX_Det(MATRIX &A)
{
	double det;

	if (A.sizeI() != A.sizeJ())	throw Common::ExceptionMemoryAllocator (FromHere(), "MATRIX_Det: Matrix should be a square matrix to calculate determinant");

	if (A.is_MatlabType() == true)
	{
		if (A.sizeI() == 1)
		{
			det = A(1, 1);
		}
		else if (A.sizeI() == 2)
		{
			det = A(1,1)*A(2,2) - A(1,2)*A(2,1);
		}
		else
		{
			det = 0.0;
#pragma omp parallel for reduction(+: det)
			for (int i = 1; i <= A.sizeI(); i++)
			{
				MATRIX A_temp = MATRIX_Minor(A, 1, i);
				double temp1 = pow(-1.0, i-1);
				double temp2 = MATRIX_Det(A_temp);

				det += temp1 * A(1, i) * temp2;
			}

		}
	}
	else
	{
		if (A.sizeI() == 1)
		{
			det = A(0, 0);
		}
		else if (A.sizeI() == 2)
		{
			det = A(0,0)*A(1,1) - A(0,1)*A(1,2);
		}
		else
		{
			det = 0.0;
#pragma omp parallel for reduction(+: det)
			for (int i = 0; i <= A.sizeI()-1; i++)
			{
				MATRIX A_temp = MATRIX_Minor(A, 0, i);
				double temp1 = pow(-1.0, i);
				double temp2 = MATRIX_Det(A_temp);

				det += temp1 * A(0, i) * temp2;
			}

		}
	}

	return (det);

}


MATRIX MATRIX_Inv(MATRIX &A)
{
	if (A.sizeI() != A.sizeJ())	throw Common::ExceptionMemoryAllocator (FromHere(), "MATRIX_Inv: Matrix should be a square matrix to calculate determinant");

	MATRIX A_inv(A.sizeI(), A.sizeI(), A.is_MatlabType());

	int 	i, j;
	int 	info;
	int 	*ipiv;



	info			= 0;
	ipiv			= new int [A.sizeI()+1];
	int 	lwork	= A_inv.sizeI() * A_inv.sizeI();

	double *work	= new double [lwork];
	double *element = new double [A.sizeI()*A.sizeI()];

	// Store element in Column-major order
	if (A.is_MatlabType() == true)
	{
		for (i = 0; i <= A.sizeI()-1; i++)
			for (j = 0; j <= A.sizeJ()-1; j++)
				element[i + j*A.sizeI()] = A(i+1,j+1);
	}
	else
	{
		for (i = 0; i <= A.sizeI()-1; i++)
			for (j = 0; j <= A.sizeJ()-1; j++)
				element[i + j*A.sizeI()] = A(i,j);
	}

	int m	= A_inv.sizeI();
	dgetrf(&m, &m, element, &m, ipiv, &info);
	dgetri(&m, element, &m, ipiv, work, &lwork, &info);

	int k;
	k	= 0;
	if (A_inv.is_MatlabType() == true)
	{
		for (j = 0; j <= A.sizeI()-1; j++)
		{
			for (i = 0; i <= A.sizeI()-1; i++)
			{
				A_inv(i+1,j+1) = element[k];
				k++;
			}
		}
	}
	else
	{
		for (j = 0; j <= A.sizeI()-1; j++)
		{
			for (i = 0; i <= A.sizeI()-1; i++)
			{
				A_inv(i,j) = element[k];
				k++;
			}
		}
	}

	delete[] ipiv;
	delete[] work;
	delete[] element;

	return (A_inv);
}


MATRIX MATRIX_Confactor(MATRIX &A)
{
	if (A.sizeI() != A.sizeJ())	throw Common::ExceptionMemoryAllocator (FromHere(), "MATRIX_Inv: Matrix should be a square matrix to calculate determinant");

	MATRIX res(A.sizeI(), A.sizeI(), A.is_MatlabType());


	if (A.is_MatlabType() == true)
	{
#pragma omp parallel
		{
			for (int i = 1; i <= A.sizeI(); i++)
			{
				for (int j = 1; j <= A.sizeI(); j++)
				{
					MATRIX ap = MATRIX_Minor(A, i, j);

					double d = MATRIX_Det(ap);
					res(i,j) = pow(-1.0, i+j) * d;
				}
			}
		}
	}
	else
	{
#pragma omp parallel
		{
			for (int i = 0; i <= A.sizeI()-1; i++)
			{
				for (int j = 0; j <= A.sizeI()-1; j++)
				{
					MATRIX ap = MATRIX_Minor(A, i, j);

					double d = MATRIX_Det(ap);
					res(i,j) = pow(-1.0, i+j) * d;
				}
			}
		}

	}

	return res;
}


MATRIX MATRIX_Adjoint(MATRIX &A)
{
	if (A.sizeI() != A.sizeJ())	throw Common::ExceptionMemoryAllocator (FromHere(), "MATRIX_Inv: Matrix should be a square matrix to calculate determinant");

	MATRIX res(A.sizeI(), A.sizeI(), A.is_MatlabType());
	MATRIX ap	= MATRIX_Confactor(A);


	if (A.is_MatlabType() == true)
	{
#pragma omp parallel
		{
			for (int i = 1; i <= A.sizeI(); i++)
			{
				for (int j = 1; j <= A.sizeI(); j++)
				{
					res(i, j) = ap(j, i);
				}
			}
		}
	}
	else
	{
#pragma omp parallel
		{
			for (int i = 0; i <= A.sizeI()-1; i++)
			{
				for (int j = 0; j <= A.sizeI()-1; j++)
				{
					res(i, j) = ap(j, i);
				}
			}
		}
	}

	return res;
}


MATRIX MATRIX_Inv2(MATRIX &A)
{
	if (A.sizeI() != A.sizeJ())	throw Common::ExceptionMemoryAllocator (FromHere(), "MATRIX_Inv: Matrix should be a square matrix to calculate determinant");
	MATRIX res(A.sizeI(), A.sizeI(), A.is_MatlabType());


	double 	d = MATRIX_Det(A);

	if (A.sizeI() == 1)		// 1 X 1 MATRIX
	{
		if (A.is_MatlabType() == true)	res(1,1) = 1/d;
		else							res(0,0) = 1/d;
	}
	else if (A.sizeI() == 2)	// 2 X 2 MATRIX
	{
		if (A.is_MatlabType() == true)
		{
			res(1,1) =  A(2,2) / d;
			res(1,2) = -A(1,2) / d;

			res(2,1) = -A(2,1) / d;
			res(2,2) =  A(1,1) / d;

		}
		else
		{
			res(0,0) =  A(1,1) / d;
			res(0,1) = -A(0,1) / d;

			res(1,0) = -A(1,0) / d;
			res(1,1) =  A(0,0) / d;
		}
	}
	else				// MORE THAN 3 X 3 MATRIX
	{
		// IT USES ADJOINT METHOD
		MATRIX ai = MATRIX_Adjoint(A);


		if (A.is_MatlabType() == true)
		{
#pragma omp parallel
			{
			for (int i = 1; i <= A.sizeI(); i++)
				for (int j = 1; j <= A.sizeI(); j++)
					res(i,j)	= ai(i,j) / d;
			}
		}
		else
		{
#pragma omp parallel
			{
			for (int i = 0; i <= A.sizeI()-1; i++)
				for (int j = 0; j <= A.sizeI()-1; j++)
					res(i,j)	= ai(i,j) / d;
			}
		}
	}

	return res;
}


}
}
