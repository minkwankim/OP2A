/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Mar 24, 2015
 *      			Author: Minkwan Kim
 *
 * CMatrix_minor.cpp
 * 			-  
 *  
 */

#include "../include/CMatrix.hpp"



//::M-07 - minor
//		 - MINOR FUNCTION
CMatrix_ver2 CMatrix_ver2::minor_mk(const int r, const int c)
{
	CMatrix_ver2	res;
	if (r >= 1 && r <= M && c >= 1 && c <= N)	res.size(M-1, N-1);		// Create new matrix smaller size
	else										Exception("[ERROR][CMatrix]: Matrix size should be bigger to perform minor_mk()");



	int  k1, k2;

	for (int i = 1; i < r; i++)
	{
		for (int j = 1; j < c; j++)
		{
			k1	= Mat_index_2D_to_1D(i, j, M-1, N-1);
			k2	= Mat_index_2D_to_1D(i, j, M, N);

			res.elements[k1]	= elements[k2];
		}

		for (int j = c+1; j <= M; j++)
		{
			k1	= Mat_index_2D_to_1D(i, j-1,	M-1, N-1);
			k2	= Mat_index_2D_to_1D(i, j,		M, N);

			res.elements[k1]	= elements[k2];
		}
	}

	for (int i = r+1; i <= M; i++)
	{
		for (int j = 1; j < c; j++)
		{
			k1	= Mat_index_2D_to_1D(i-1,	j, M-1, N-1);
			k2	= Mat_index_2D_to_1D(i, 	j, M, N);

			res.elements[k1]	= elements[k2];
		}

		for (int j = c+1; j <= M; j++)
		{
			k1	= Mat_index_2D_to_1D(i-1,	j-1,	M-1, N-1);
			k2	= Mat_index_2D_to_1D(i, 	j,		M, N);

			res.elements[k1]	= elements[k2];
		}
	}

	return res;
}

void CMatrix_ver2::minor_mk2(const int r, const int c)
{
	delete_row(r);
	delete_column(c);
}
