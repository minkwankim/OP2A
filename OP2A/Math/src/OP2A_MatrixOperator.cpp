/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 26, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_MatrixOperator.cpp
 * 			-  
 *  
 */




#include "Math/include/OP2A_Matrix.hpp"
#include "Common/include/Exception_MemoryAllocator.hpp"
#include "Common/include/Exception_DimensionMatch.hpp"


namespace OP2A{
namespace Math{



double& MATRIX::operator() (const unsigned int i, const unsigned int j)
{
	if (m_matlabType == true)
	{
		if (i > m_I)	throw Common::ExceptionMemoryAllocator(FromHere(), "Exceed the row size of a MATRIX");
		if (i == 0) 	throw Common::ExceptionMemoryAllocator(FromHere(), "Matlab-Type MATRIX Index should be non-zero integer");

		if (j > m_J)	throw Common::ExceptionMemoryAllocator(FromHere(), "Exceed the row size of a MATRIX");
		if (j == 0) 	throw Common::ExceptionMemoryAllocator(FromHere(), "Matlab-Type MATRIX Index should be non-zero integer");

		return (m_data[i-1][j-1]);
	}
	else
	{
		if (i >= m_I)	throw Common::ExceptionMemoryAllocator(FromHere(), "Exceed the row size of a C-type MATRIX");
		if (j >= m_J)	throw Common::ExceptionMemoryAllocator(FromHere(), "Exceed the row size of a C-type MATRIX");

		return (m_data[i][j]);
	}
}


void MATRIX::operator*= (const double s)
{
#pragma omp parallel for
	for (int i = 0; i <= m_I-1; i++)
	{
#pragma ivdep
		for (int j = 0; j <= m_J-1; j++)
		{
			m_data[i][j] *= s;
		}
	}
}

void MATRIX::operator/= (const double s)
{
#pragma omp parallel for
	for (int i = 0; i <= m_I-1; i++)
	{
#pragma ivdep
		for (int j = 0; j <= m_J-1; j++)
		{
			m_data[i][j] /= s;
		}
	}
}

void MATRIX::operator+= (const double s)
{
#pragma omp parallel for
	for (int i = 0; i <= m_I-1; i++)
	{
#pragma ivdep
		for (int j = 0; j <= m_J-1; j++)
		{
			m_data[i][j] += s;
		}
	}
}

void MATRIX::operator-= (const double s)
{
#pragma omp parallel for
	for (int i = 0; i <= m_I-1; i++)
	{
#pragma ivdep
		for (int j = 0; j <= m_J-1; j++)
		{
			m_data[i][j] -= s;
		}
	}
}


void MATRIX::operator+= (const MATRIX & S)
{
	if (m_I != S.m_I)	throw Common::ExceptionMemoryAllocator(FromHere(), "Operator +=: Two matrix should have same dimension");
	if (m_J != S.m_J)	throw Common::ExceptionMemoryAllocator(FromHere(), "Operator +=: Two matrix should have same dimension");


#pragma omp parallel for
	for (int i = 0; i <= m_I-1; i++)
	{
#pragma ivdep
		for (int j = 0; j <= m_J-1; j++)
		{
			m_data[i][j] += S.m_data[i][j];
		}
	}
}

void MATRIX::operator-= (const MATRIX & S)
{
	if (m_I != S.m_I)	throw Common::ExceptionMemoryAllocator(FromHere(), "Operator -=: Two matrix should have same dimension");
	if (m_J != S.m_J)	throw Common::ExceptionMemoryAllocator(FromHere(), "Operator -=: Two matrix should have same dimension");


#pragma omp parallel for
	for (int i = 0; i <= m_I-1; i++)
	{
#pragma ivdep
		for (int j = 0; j <= m_J-1; j++)
		{
			m_data[i][j] -= S.m_data[i][j];
		}
	}
}




MATRIX& MATRIX::operator= (std::vector< std::vector<double> >& A)
{
	if (m_I !=  A.size() || m_J != A[0].size())	resize(A.size(), A[0].size());

	m_data	= A;
	return *this;
}



/*
 * Operator Functions
 */
MATRIX operator- (MATRIX &A, MATRIX &B)
{
	if (A.sizeI() != B.sizeJ() || A.sizeJ() != B.sizeJ())	throw Common::ExceptionDimensionMatch (FromHere(), "MATRIX Operator - :Two matrix should have same dimension");
	if (A.is_MatlabType() != B.is_MatlabType())				throw Common::ExceptionDimensionMatch (FromHere(), "MATRIX Operator - :Two matrix should have same Type");

	MATRIX C(A.sizeI(), A.sizeJ(), A.is_MatlabType());

	if (C.is_MatlabType() == true)
	{
#pragma omp parallel for
		for (int i = 1; i <= C.sizeI(); i++)
		{
#pragma ivdep
			for (int j = 1; j <= C.sizeJ(); j++)
			{
				C(i,j) = A(i,j) - B(i,j);
			}
		}
	}
	else
	{
#pragma omp parallel for
		for (int i = 0; i <= C.sizeI()-1; i++)
		{
#pragma ivdep
			for (int j = 0; j <= C.sizeJ()-1; j++)
			{
				C(i,j) = A(i,j) - B(i,j);
			}
		}
	}

	return (C);
}


MATRIX operator+ (MATRIX &A, MATRIX &B)
{
	if (A.sizeI() != B.sizeJ() || A.sizeJ() != B.sizeJ())	throw Common::ExceptionDimensionMatch (FromHere(), "MATRIX Operator - :Two matrix should have same dimension");
	if (A.is_MatlabType() != B.is_MatlabType())				throw Common::ExceptionDimensionMatch (FromHere(), "MATRIX Operator - :Two matrix should have same Type");

	MATRIX C(A.sizeI(), A.sizeJ(), A.is_MatlabType());

	if (C.is_MatlabType() == true)
	{
#pragma omp parallel for
		for (int i = 1; i <= C.sizeI(); i++)
		{
#pragma ivdep
			for (int j = 1; j <= C.sizeJ(); j++)
			{
				C(i,j) = A(i,j) + B(i,j);
			}
		}
	}
	else
	{
#pragma omp parallel for
		for (int i = 0; i <= C.sizeI()-1; i++)
		{
#pragma ivdep
			for (int j = 0; j <= C.sizeJ()-1; j++)
			{
				C(i,j) = A(i,j) + B(i,j);
			}
		}
	}

	return (C);
}



MATRIX operator* (double a, MATRIX &A)
{
	MATRIX C(A.sizeI(), A.sizeJ(), A.is_MatlabType());

	if (C.is_MatlabType() == true)
	{
#pragma omp parallel for
		for (int i = 1; i <= C.sizeI(); i++)
		{
#pragma ivdep
			for (int j = 1; j <= C.sizeJ(); j++)
			{
				C(i,j) = a * A(i,j);
			}
		}
	}
	else
	{
#pragma omp parallel for
		for (int i = 0; i <= C.sizeI()-1; i++)
		{
#pragma ivdep
			for (int j = 0; j <= C.sizeJ()-1; j++)
			{
				C(i,j) = a * A(i,j);
			}
		}
	}

	return (C);
}


MATRIX operator* (MATRIX &A, MATRIX &B)
{
	if (A.sizeJ() != B.sizeI())	throw Common::ExceptionDimensionMatch (FromHere(), "MATRIX Operator * : Check Dimension of matrix");
	if (A.is_MatlabType() != B.is_MatlabType())				throw Common::ExceptionDimensionMatch (FromHere(), "MATRIX Operator - :Two matrix should have same Type");

	MATRIX C(A.sizeI(), B.sizeJ(), A.is_MatlabType());


	if (C.is_MatlabType() == true)
	{
#pragma ivdep
		for (int i = 1; i <= C.sizeI(); i++)
		{
#pragma ivdep
			for (int j = 1; j <= C.sizeJ(); j++)
			{

				double sum_temp = 0.0;
//#pragma omp parallel for reduction(+: sum_temp)
				for (int index = 1; index <= A.sizeJ(); index++)
				{
					sum_temp += A(i, index) * B(index, j);
				}

				C(i, j)	= sum_temp;
			}
		}
	}
	else
	{
#pragma ivdep
		for (int i = 0; i <= C.sizeI()-1; i++)
		{

#pragma ivdep
			for (int j = 0; j <= C.sizeJ()-1; j++)
			{
				double sum_temp = 0.0;
//#pragma omp parallel for reduction(+: sum_temp)
				for (int index = 0; index <= A.sizeJ()-1; index++)
				{
					sum_temp += A(i, index) * B(index, j);
				}

				C(i, j)	= sum_temp;
			}
		}
	}

	return (C);
}


std::vector<double> operator* (MATRIX &A, std::vector<double> &B)
{
	if (A.sizeJ() != B.size())	throw Common::ExceptionDimensionMatch (FromHere(), "MATRIX Operator * : Check Dimension of matrix and vector");

	std::vector<double>	res(A.sizeI(), 0.0);

	if (A.is_MatlabType() == true)
	{
#pragma omp parallel for
		for (int i = 1; i <= A.sizeI(); i++)
		{
			double temp = 0.0;
//#pragma omp parallel for reduction(+: temp)
			for (int j = 1; j <= B.size(); j++)
			{
				temp += A(i,j) * B[j-1];
			}

			res[i-1] = temp;
		}
	}
	else
	{
#pragma omp parallel for
		for (int i = 0; i <= A.sizeI()-1; i++)
		{
			double temp = 0.0;
//#pragma omp parallel for reduction(+: temp)
			for (int j = 0; j <= B.size()-1; j++)
			{
				temp += A(i,j) * B[j];
			}

			res[i] = temp;
		}
	}

	return (res);
}




}
}
