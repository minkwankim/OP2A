/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Mar 24, 2015
 *      			Author: Minkwan Kim
 *
 * CMatrix_INV.cpp
 * 			-  
 *  
 */




#include "../include/CMatrix.hpp"
CMatrix_ver2	CMatrix_INV(CMatrix_ver2& A)
{
	int 	i, j;
	int 	info;
	int 	*ipiv;
	int		m;
	CMatrix_ver2 A_inv;

	if (A.M != A.N)	Exception("[ERROR][matrix_MK]: Square matrix is required to calculate inverse matrix");

	m	= A.M;
	A_inv.size(m, m);

	info			= 0;
	ipiv			= new int [A.M+1];
	int 	lwork	= A_inv.M*A_inv.M;
	double *work	= new double [lwork];
	double *element = new double [m*m];

	// Store element in Column-major order
	for (i = 0; i <= m-1; i++)
		for (j = 0; j <= m-1; j++)
			element[i + j*m] = A(i+1,j+1);

	dgetrf(&A_inv.M, &A_inv.N, element, &A_inv.M, ipiv, &info);
	dgetri(&A_inv.M, element, &A_inv.M, ipiv, work, &lwork, &info);

	int k;
	k	= 0;
	for (j = 0; j <= m-1; j++)
	{
		for (i = 0; i <= m-1; i++)
		{
			A_inv(i+1,j+1) = element[k];
			k++;
		}
	}

	delete[] ipiv;
	delete[] work;
	delete[] element;

	return (A_inv);
}




// CF-04:: Determinant
double	CMatrix_DET(CMatrix_ver2& A)
{
	int row, col;
	int c;
	double temp1, temp2;
	double d;			// value of the determinant

	d = 0;
	row = A.M;
	col = A.N;


	if (row != col)
	{
		throw Exception("[ERROR][matrix_MK]:  Error using ==> det.  matrix must be square.");
	}
	else
	{
		if (row == 1)
	    {
			d = A(1, 1);
		}
		else if (row == 2)
		{
			d = A(1,1) * A(2,2) - A(1,2) * A(2,1);
		}
		else
		{
			for (c = 1; c <= col; c++)
			{
				CMatrix_ver2 M	= A.minor_mk(1, c);
				temp1 			= pow(-1.0, c);
				temp2 			= CMatrix_DET(M);

				d += temp1 * A(1, c) * temp2;  // faster than with pow()
			}
		}
	}

	return d;
}
