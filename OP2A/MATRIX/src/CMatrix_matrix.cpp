/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Mar 24, 2015
 *      			Author: Minkwan Kim
 *
 * CMatrix_matrix.cpp
 * 			-  
 *  
 */


#include "../include/CMatrix.hpp"


void CMatrix_ver2::diag()
{
	if (N*M == 0)	throw Exception("[ERROR][CMatrix]-[diag()]: Size is not defined yet.");
	if (M != N)		throw Exception("[ERROR][CMatrix]-[diag()]: Diagonal natrix function is only applicable for a square matrix.");

	int k;
	for (int i = 1; i <= M; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			k	= Mat_index_2D_to_1D(i, j, M, N);
			if (i == j)	elements[k] = 1.0;
			else		elements[k] = 0.0;
		}
	}
}

void CMatrix_ver2::diag(int n)
{
	if (N*M == 0)	size(n, n);

	if (M != n || N != n)
	{
		delete [] elements;
		M = n;
		N = n;
		elements	= new double [M*N];
	}

	int k;
	for (int i = 1; i <= M; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			k	= Mat_index_2D_to_1D(i, j, M, N);
			if (i == j)	elements[k] = 1.0;
			else		elements[k] = 0.0;
		}
	}
}


void CMatrix_ver2::ones()
{
	if (N*M == 0)	throw Exception("[ERROR][CMatrix]-[ones()]: Size is not defined yet.");
	for (int i = 0; i <= M*N-1; i++)	elements[i]	= 1.0;
}

void CMatrix_ver2::ones(int m, int n)
{
	if (N*M == 0)	size(m, n);

	if (M != m || N != n)
	{
		delete [] elements;
		M = m;
		N = n;
		elements	= new double [M*N];
	}

	for (int i = 0; i <= M*N-1; i++)	elements[i]	= 1.0;
}


void CMatrix_ver2::zeros()
{
	if (N*M == 0)	throw Exception("[ERROR][CMatrix]-[zeros()]: Size is not defined yet.");
	for (int i = 0; i <= M*N-1; i++)	elements[i]	= 0.0;
}

void CMatrix_ver2::zeros(int m, int n)
{
	if (N*M == 0)	size(m, n);

	if (M != m || N != n)
	{
		delete [] elements;
		M = m;
		N = n;
		elements	= new double [M*N];
	}

	for (int i = 0; i <= M*N-1; i++)	elements[i]	= 0.0;
}



CMatrix_ver2	CMatrix_ZEROS(int i, int j)
{
	CMatrix_ver2 A;
	A.zeros(i, j);
	return A;
}

CMatrix_ver2 	CMatrix_ONES(int i, int j)
{
	CMatrix_ver2 A;
	A.ones(i, j);
	return A;
}

CMatrix_ver2 	CMatrix_DIAG(int i)
{
	CMatrix_ver2 A;
	A.diag(i);

	return A;
}



