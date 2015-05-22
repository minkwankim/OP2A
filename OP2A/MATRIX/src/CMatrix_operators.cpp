/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Mar 23, 2015
 *      			Author: Minkwan Kim
 *
 * CMatrix_operators.cpp
 * 			-  
 *  
 */


#include "../include/CMatrix.hpp"



double& CMatrix_ver2::operator()(const int i, const int j)
{
	if (i <= 0 || i > M)	throw Exception("[ERROR][CMatrix]:  i-Subscript indices must either be real positive integers or logicals and it should be with a matrix size.");
	if (j <= 0 || j > N)	throw Exception("[ERROR][CMatrix]:  j-Subscript indices must either be real positive integers or logicals and it should be with a matrix size.");

	int k	= Mat_index_2D_to_1D(i, j, M, N);
	return (elements[k]);
}


CMatrix_ver2& CMatrix_ver2::operator= (const CMatrix_ver2 &A)
{
	int i;

	if (N*M == 0)
	{
		size(A.M, A.N);
		for (int i = 0; i <= N*M-1; i++)	elements[i] = A.elements[i];
		return *this;
	}

	if (M != A.M || N != A.N)	resize(A.M, A.N);

	for (int i = 0; i <= N*M-1; i++)	elements[i] = A.elements[i];
	return *this;
}

CMatrix_ver2& CMatrix_ver2::operator= (const vector <vector <double> > &A)
{
	int m	= A.size();
	int n	= A[0].size();
	int k;


	if (N*M == 0)
	{
		size(m, n);

		k	= 0;
		for (int i = 0; i <= M-1; i++)
		{
			for (int j = 0; j <= N-1; j++)
			{
				elements[k]	= A[i][j];
				k++;
			}
		}

		return *this;
	}

	k	= 0;
	for (int i = 0; i <= M-1; i++)
	{
		for (int j = 0; j <= N-1; j++)
		{
			elements[k]	= A[i][j];
			k++;
		}
	}

	return *this;
}


CMatrix_ver2& CMatrix_ver2::operator+= (const double a)
{
	if (N*M == 0)	throw Exception("[ERROR][CMatrix]:  Matrix size is not assigned yet.");
	for (int i = 0; i <= N*M-1; i++)	elements[i] += a;
	return *this;
}


CMatrix_ver2& CMatrix_ver2::operator+= (const CMatrix_ver2& A)
{
	if (N*M == 0)				throw Exception("[ERROR][CMatrix]:  Matrix size is not assigned yet.");
	if (M != A.M || N != A.N)	throw Exception("[ERROR][CMatrix]:  Two matrix should have same size.");

	for (int i = 0; i <= N*M-1; i++)	elements[i] += A.elements[i];
	return *this;
}


CMatrix_ver2& CMatrix_ver2::operator-= (const double a)
{
	if (N*M == 0)	throw Exception("[ERROR][CMatrix]:  Matrix size is not assigned yet.");
	for (int i = 0; i <= N*M-1; i++)	elements[i] -= a;
	return *this;
}


CMatrix_ver2& CMatrix_ver2::operator-= (const CMatrix_ver2& A)
{
	if (N*M == 0)				throw Exception("[ERROR][CMatrix]:  Matrix size is not assigned yet.");
	if (M != A.M || N != A.N)	throw Exception("[ERROR][CMatrix]:  Two matrix should have same size.");

	for (int i = 0; i <= N*M-1; i++)	elements[i] -= A.elements[i];
	return *this;
}







CMatrix_ver2	operator+ (const CMatrix_ver2& A, 	const CMatrix_ver2& B)
{
	CMatrix_ver2 result;

	if (A.M == A.M && A.N == B.N)
	{
		result.size(A.M, B.N);
		for (int i = 0; i <= A.M*A.N-1; i++)	result.elements[i] = A.elements[i] + B.elements[i];
	}
	else
	{
		Exception("[ERROR][matrix_MK]: Dimensions do not match.");
	}

	return result;
}

CMatrix_ver2	operator+ (const CMatrix_ver2& A, 	const double B)
{
	CMatrix_ver2 result	= A;
	result	+= B;

	return result;
}

CMatrix_ver2	operator+ (const double B, 			const CMatrix_ver2& A)
{
	CMatrix_ver2 result	= A;
	result	+= B;

	return result;
}



CMatrix_ver2	operator- (const CMatrix_ver2& A, 	const CMatrix_ver2& B)
{
	CMatrix_ver2 result;

	if (A.M == A.M && A.N == B.N)
	{
		result.size(A.M, B.N);
		for (int i = 0; i <= A.M*A.N-1; i++)	result.elements[i] = A.elements[i] - B.elements[i];
	}
	else
	{
		Exception("[ERROR][matrix_MK]: Dimensions do not match.");
	}

	return result;
}

CMatrix_ver2	operator- (const CMatrix_ver2& A, 	const double B)
{
	CMatrix_ver2 result	= A;

	for (int i = 0; i <= A.M*A.N-1; i++)	result.elements[i] = A.elements[i] - B;
	return result;
}

CMatrix_ver2	operator- (const double B, 			const CMatrix_ver2& A)
{
	CMatrix_ver2 result	= A;

	for (int i = 0; i <= A.M*A.N-1; i++)	result.elements[i] = B - A.elements[i];
	return result;
}

CMatrix_ver2	operator- (const CMatrix_ver2& A)
{
	CMatrix_ver2 result	= A;
	for (int i = 0; i <= A.M*A.N-1; i++)	result.elements[i] = -A.elements[i];
	return result;
}


CMatrix_ver2	operator* (const CMatrix_ver2& A, 	const CMatrix_ver2& B)
{
	int m, n, k;
	double alpha, beta;

	CMatrix_ver2 res;
	alpha	= 1.0;
	beta	= 0.0;
	if (A.N == B.M)
	{
		m	= A.M;
		k	= A.N;
		n	= B.N;

		res.size(m, n);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A.elements, k, B.elements, n, beta, res.elements, n);
	}
	else
	{
		Exception("[ERROR][matrix_MK]: Dimensions do not match.");
	}
	return res;
}


CMatrix_ver2	operator* (const double B,			const CMatrix_ver2& A)
{
	CMatrix_ver2 res = A;
	res.mul(B);
	return res;
}
