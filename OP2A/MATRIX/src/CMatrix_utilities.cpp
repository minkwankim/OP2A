/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Mar 24, 2015
 *      			Author: Minkwan Kim
 *
 * CMatrix_utilities.cpp
 * 			-  
 *  
 */




#include "../include/CMatrix.hpp"


void CMatrix_ver2::resize(int m, int n)
{
	delete [] elements;
	M = m;
	N = n;
	elements	= new double [M*N];
}


void CMatrix_ver2::add_column()
{
	double *temp	= new double [M*(N+1)];

	int k1;
	int k2;

	for (int i = 1; i <= M; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			k1	= Mat_index_2D_to_1D(i, j, M, N+1);
			k2	= Mat_index_2D_to_1D(i, j, M, N);

			temp[k1]	= elements[k2];
		}

		k1	= Mat_index_2D_to_1D(i, N+1, M, N+1);
		temp[k1]	= 0.0;
	}

	delete [] elements;

	N = N + 1;
	elements	= temp;
}


void CMatrix_ver2::add_column(int a)
{
	double *temp	= new double [M*(N+1)];

	int k1;
	int k2;

	for (int i = 1; i <= M; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			if (j < a)	k1	= Mat_index_2D_to_1D(i, j,		M, N+1);
			else		k1	= Mat_index_2D_to_1D(i, j+1, 	M, N+1);

			k2	= Mat_index_2D_to_1D(i, j, M, N);

			temp[k1]	= elements[k2];
		}

		k1	= Mat_index_2D_to_1D(i, a, M, N+1);
		temp[k1]	= 0.0;
	}

	delete [] elements;

	N = N + 1;
	elements	= temp;
}


void CMatrix_ver2::add_row()
{
	double *temp	= new double [(M+1)*N];

	int k1;
	int k2;

	for (int i = 1; i <= M; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			k1	= Mat_index_2D_to_1D(i, j, M+1, N);
			k2	= Mat_index_2D_to_1D(i, j, M, 	N);

			temp[k1]	= elements[k2];
		}
	}

	for (int j = 1; j <= N; j++)
	{
		k1	= Mat_index_2D_to_1D(N+1, j, M+1, N);
		temp[k1]	= 0.0;
	}

	delete [] elements;

	M = M + 1;
	elements	= temp;
}



void CMatrix_ver2::add_row(int a)
{
	double *temp	= new double [(M+1)*N];

	int k1;
	int k2;

	for (int i = 1; i <= M; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			if (i < a)	k1	= Mat_index_2D_to_1D(i, 	j,	M+1, N);
			else		k1	= Mat_index_2D_to_1D(i+1, 	j, 	M+1, N);

			k2	= Mat_index_2D_to_1D(i, j, M, N);

			temp[k1]	= elements[k2];
		}
	}

	for (int j = 1; j <= N; j++)
	{
		k1	= Mat_index_2D_to_1D(a, j, M+1, N);
		temp[k1]	= 0.0;
	}



	delete [] elements;

	M = M + 1;
	elements	= temp;
}




void CMatrix_ver2::delete_column(int a)
{
	double *temp	= new double [M*(N-1)];

	int k1;
	int k2;

	for (int j = 1; j < a; j++)
	{
		for (int i = 1; i <= M; i++)
		{
			k1	= Mat_index_2D_to_1D(i, j, M, N-1);
			k2	= Mat_index_2D_to_1D(i, j, M, N);
			temp[k1]	= elements[k2];
		}
	}

	for (int j = a; j <= N-1; j++)
	{
		for (int i = 1; i <= M; i++)
		{
			k1	= Mat_index_2D_to_1D(i, j, 	 M, N-1);
			k2	= Mat_index_2D_to_1D(i, j+1, M, N);
			temp[k1]	= elements[k2];
		}
	}

	delete [] elements;
	N = N - 1;
	elements	= temp;
}



void CMatrix_ver2::delete_row(int a)
{
	double *temp	= new double [(M-1)*N];

	int k1;
	int k2;

	for (int i = 1; i < a; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			k1	= Mat_index_2D_to_1D(i, j, M-1, N);
			k2	= Mat_index_2D_to_1D(i, j, M, 	N);
			temp[k1]	= elements[k2];
		}
	}

	for (int i = a; i <= M-1; i++)
	{
		for (int j = 1; j <= N; j++)
		{
			k1	= Mat_index_2D_to_1D(i,		j,	M-1, 	N);
			k2	= Mat_index_2D_to_1D(i+1,	j,	M, 		N);
			temp[k1]	= elements[k2];
		}
	}

	delete [] elements;
	M = M - 1;
	elements	= temp;
}
