/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Mar 13, 2015
 *      			Author: Minkwan Kim
 *
 * OPPA_Matrix.cpp
 * 			-  
 *  
 */
#include "omp.h"
#include <mkl.h>
#include "../include/OPPA_Matrix.hpp"
#include "../include/matrix.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"

/*
 * ======================================
 * Constructors
 */

Matrix::Matrix()
{
	N	= 0;
	M	= 0;
}

Matrix::Matrix(int m, int n)
{
	M	= m;
	N	= n;
	elements	= vector_2D(M, N, 0.0);
}

Matrix::Matrix(const Matrix &other)
{
	M			= other.M;
	N			= other.N;
	elements	= other.elements;
}

Matrix::Matrix(vector< vector <double> > &other)
{
	M	= elements.size();
	N	= elements[0].size();
	elements	= other;
}

Matrix::~Matrix()
{

}



/*
 * =========================================
 * Internal functions
 * =========================================
 */
// Member functions
//
//::M-01 - size()
//		 - ASSIGNING THE SIZE
void Matrix::size(int m, int n)
{
	M			= m;
	N 			= n;
	elements.resize(m);
	for (int i = 0; i <= m; i++)	elements[i].resize(n);
}


//::M-02 - diag()
// 		 - DIAGONAL MATRIX
void Matrix::diag()
{
	if (N*M == 0)	throw Exception("[ERROR][Matrix_MK]-[diag()]: Size is not defined yet. or it is not a matrix");
	if (M != N)		throw Exception("[ERROR][Mmatrix_MK]-[diag()]: Diagonal natrix function is only applicable for a square matrix.");

	for (int i= 0; i <= M-1; i++)
		for (int j = 0; j <= N-1; j++)	elements[i][j]	= 0.0;

	for (int i= 0; i <= M-1; i++)	elements[i][i]	= 1.0;
}

void Matrix::diag(int n)
{
	if (N*M == 0)			size(n, n);
	if (N != n || M != n)	size(n, n);

	for (int i= 0; i <= M-1; i++)
		for (int j = 0; j <= N-1; j++)	elements[i][j]	= 0.0;

	for (int i= 0; i <= M-1; i++)	elements[i][i]	= 1.0;
}


//::M-03 - ones()
// 		 - ONES FUNCTIONS
void Matrix::ones()
{
	for (int i= 0; i <= M-1; i++)
		for (int j = 0; j <= N-1; j++)	elements[i][j]	= 1.0;
}

void Matrix::ones(int i, int j)
{
	size(i, j);

	for (int i= 0; i <= M-1; i++)
		for (int j = 0; j <= N-1; j++)	elements[i][j]	= 0.0;
}



//::M-04 - zeros()
// 		 - ZEROSS FUNCTIONS
void Matrix::zeros()
{
	for (int i= 0; i <= M-1; i++)
		for (int j = 0; j <= N-1; j++)	elements[i][j]	= 0.0;
}

void Matrix::zeros(int i, int j)
{
	size(i, j);

	for (int i= 0; i <= M-1; i++)
		for (int j = 0; j <= N-1; j++)	elements[i][j]	= 0.0;
}


//::M-05 - OPERATORS
double& Matrix::operator()(const int i, const int j)
{
	if (i < 0 || i > M)	throw Exception("[ERROR][matrix_MK]:  Subscript indices must either be real positive integers or logicals and it should be with a matrix size.");
	if (j < 0 || j > N)	throw Exception("[ERROR][matrix_MK]:  Subscript indices must either be real positive integers or logicals and it should be with a matrix size.");

	return (elements[i][j]);
}

Matrix& Matrix::operator= (const Matrix& a)
{
	int i, j;
	int flag;

	if (N*M == 0)
	{
		size(a.M, a.N);
		flag = 0;
	}
	else
	{
		flag = 1;
		if (M == a.M && N == a.N)
		{
			flag = 0;
		}
	}

	if (flag == 0)
	{
		for (i = 0; i <= M-1; i++)
			for (j = 0; j <= N-1; j++)
				elements[i][j] = a.elements[i][j];
	}
	else
	{
		Exception("[ERROR][matrix_MK]: Two matrix should have same size.");
	}

	return *this;
}


Matrix& Matrix::operator= (vector< vector<double> >& a)
{
	int flag;

	if (N*M == 0)	size(a.size(), a[0].size());

	if (M != a.size() || N != a[0].size())	size(a.size(), a[0].size());

	elements	= a;
	return *this;
}


//::M-06 - ADD / SUBRACT / MULTIPLY
//		 - add
Matrix& Matrix::add(double v)
{
	for (int i = 0; i <= M-1; i++)
		for (int j = 0; j <= N-1; j++)
			elements[i][j] = elements[i][j] + v;

	return *this;
}

Matrix& Matrix::sub(double v)
{
	for (int i = 0; i <= M-1; i++)
		for (int j = 0; j <= N-1; j++)
			elements[i][j] = elements[i][j] - v;

	return *this;
}

Matrix& Matrix::mul(double v)
{
	for (int i = 0; i <= M-1; i++)
		for (int j = 0; j <= N-1; j++)
			elements[i][j] = elements[i][j] * v;

	return *this;
}


//::M-07 - minor
//		 - MINOR FUNCTION
Matrix Matrix::minor_mk(const int r, const int c)
{
	int i, j;
	int res_r, res_c;

	Matrix	res;

	if (r >= 0 && r <= M-1 && c >= 0 && c <= N-1)
	{
		res.size(M-1, N-1);		// Create new matrix smaller size

		res_r = 0;
		for (i = 0; i <= M-1; i++)
		{
			if(i != r)
			{
				res_c = 0;
				for (j = 0; j <= N-1; j++)
				{
					if(j != c)
					{
						res.elements[res_r][res_c]	= elements[i][j];
						res_c++;
					}
				}
				res_r++;
			}
		}
	}

	return res;
}


/*
 * ===================================================================================================
 * 		Function for Matrix calcualtions
 * 												- ver 1.0
 *
 * 		[Original verion written]
 * 									- by Minkwan Kim
 * 									- on 17/Mar/2015
 * 		[last modified]
 * 									- by Minkwan Kim
 * 									- on 17/MAr/2015
 * ========================================================
 */
void check_matrix_2D_square(vector<vector<double> > &A)
{
	int M	= A.size();
	int N	= A[0].size();

	if (N != M)
	{
		Error_message_type	error;
		error.message		= "IT IS NOT A SQUARE MATRIX. PLEASE CHECK THE DIMENSION OF METRIX";
		error.module_name	= "MATRIX: OPPA_MATRIX.cpp";
		error.print_message();
	}
}

vector<vector<double> > matrix_minor(vector<vector<double> > &A, int r, int c)
{
	int M = A.size();
	int N = A[0].size();
	vector< vector<double> > res	= vector_2D(M-1, N-1, 0.0);

	if (r < 0 || r > M || c < 0 || c > N)
	{
		Error_message_type	error;
		error.message		= "Problem in the number of row and column. PLEASE CHECK THE DIMENSION OF METRIX";
		error.module_name	= "MATRIX: OPPA_MATRIX.cpp";
		error.print_message();
	}

	int res_r, res_c;

	res_r = 0;
	for (int i = 0; i <= M-1; i++)
	{
		if(i != r)
		{
			res_c = 0;
			for (int j = 0; j <= N-1; j++)
			{
				if(j != c)
				{
					res[res_r][res_c]	= A[i][j];
					res_c++;
				}
			}
			res_r++;
		}
	}

	return res;
}


double matrix_det(vector<vector<double> > &A)
{
	double d;
	check_matrix_2D_square(A);

	int N = A.size();

	if (N == 1)
	{
		d = A[0][0];
	}
	else if (N == 2)
	{
		d = A[0][0]*A[1][1] - A[0][1]*A[1][0];
	}
	else
	{
		for (int c = 0; c <= N-1; c++)
		{
			vector< vector<double> > M	=  matrix_minor(A, 0, c);
			double temp1 = pow(-1.0, c);
			double temp2 = matrix_det(M);

			d += temp1 * A[0][c] * temp2;  // faster than with pow()
		}
	}

	return d;
}



vector<vector<double> > matrix_confactor(vector <vector<double> > &A)
{
	check_matrix_2D_square(A);
	int M = A.size();

	vector <vector <double> > res	= vector_2D(M, M, 0.0);

	for (int i = 0; i <= M-1; i++)
	{
		for (int j = 0; j <= M-1; j++)
		{
			vector <vector <double> > ap = matrix_minor(A, i, j);
			double d = matrix_det(ap);
			res[i][j] = pow(-1.0, i+j) * d;
		}
	}

	return res;
}


vector<vector<double> > matrix_adjoint(vector<vector<double> >& A)
{
	int M	= A.size();
	vector <vector <double> > res	= vector_2D(M, M, 0.0);
	vector <vector <double> > ap 	= matrix_confactor(A);


	for (int i = 0; i <= M-1; i++)
	{
		for (int j = 0; j <= M-1; j++)
		{
			res[i][j] = ap[j][i];
		}
	}

	return res;
}


vector<vector<double> > matrix_inv(vector<vector<double> > &A)
{
	check_matrix_2D_square(A);

	int 	M = A.size();
	vector< vector<double> > res	= vector_2D(M, M, 0.0);

	double 	d = matrix_det(A);

	if (M == 1)		// 1 X 1 MATRIX
	{
		res[0][0] = 1/d;
	}
	else if (M == 2)	// 2 X 2 MATRIX
	{
		res[0][0] =  A[1][1] / d;
		res[0][1] = -A[0][1] / d;

		res[1][0] = -A[1][0] / d;
		res[1][1] =  A[0][0] / d;
	}
	else				// MORE THAN 3 X 3 MATRIX
	{
		// IT USES ADJOINT METHOD
		vector <vector <double> > ai = matrix_adjoint(A);

		for (int i = 0; i <= M-1; i++)
			for (int j = 0; j <= M-1; j++)
				res[i][j]	= ai[i][j] / d;
	}

	return res;
}

vector<vector<double> > matrix_inv2(vector<vector<double> > &A)
{
	int 	M = A.size();

	vector< vector<double> > res	= vector_2D(M, M, 0.0);
	CMatrix	A_temp;	A_temp.size(M, M);

#pragma omp parallel for num_threads(M)
	for (int i = 0; i <= M-1; i++)
		for (int j = 0; j <= M-1; j++)
			A_temp(i+1, j+1)	= A[i][j];

	A_temp	= INV(A_temp);

#pragma omp parallel for num_threads(M)
	for (int i = 0; i <= M-1; i++)
		for (int j = 0; j <= M-1; j++)
			res[i][j]	= A_temp(i+1, j+1);

	return res;
}
















vector<double> matrix_multi_mat_vec(vector<vector<double> > &A, vector<double> &X)
{
	int N	= A.size();
	int M 	= A[0].size();

	if (N != X.size())
	{
		Error_message_type	error;
		error.message		= "Problem in the size of Matrix for multiplification. PLEASE CHECK THE DIMENSION OF METRIX";
		error.module_name	= "MATRIX: OPPA_MATRIX.cpp";
		error.print_message();
	}

	vector<double> B(M, 0.0);

	for (int i = 0; i <= M-1; i++)
	{
		double aux = 0.0;
		for (int j = 0; j <= M-1; j++)	aux	+= A[i][j]*X[j];

		B[i]	= aux;
	}

	return(B);
}



vector< vector<double> > matrix_multi_mat_mat(vector<vector<double> > &A, vector<vector<double> > &B, int nt)
{
	kmp_set_defaults("KMP_AFFINITY = scatter");

	int M	= A.size();
	int N 	= A[0].size();
	int L	= B[0].size();

	if (N != B.size())
	{
		Error_message_type	error;
		error.message		= "Problem in the size of Matrix for multiplification. PLEASE CHECK THE DIMENSION OF METRIX";
		error.module_name	= "MATRIX: OPPA_MATRIX.cpp";
		error.print_message();
	}


	vector<vector<double> > C	= vector_2D(M, L, 0.0);

	for (int i = 0; i <= M-1; i++)
	{
#pragma omp parallel for num_threads(nt)
		for (int j = 0; j <= L-1; j++)
		{
			double aux = 0.0;
			for (int k = 0; k <= N-1; k++)	aux	+= A[i][k]*B[k][j];

			C[i][j]	= aux;
		}
	}

	return(C);
}



vector< vector<double> > matrix_sub_mat_mat(vector<vector<double> > &A, vector<vector<double> > &B, int nt)
{
	kmp_set_defaults("KMP_AFFINITY = scatter");

	int M	= A.size();
	int N 	= A[0].size();

	vector<vector<double> > C	= vector_2D(M, N, 0.0);

	for (int i = 0; i <= M-1; i++)
	{
#pragma omp parallel for num_threads(nt)
		for (int j = 0; j <= N-1; j++)
		{
			C[i][j]	= A[i][j]	- B[i][j];
		}
	}

	return(C);

}

vector<double> matrix_sub_vec_vec(vector<double> &A, vector<double> &B, int nt)
{
	kmp_set_defaults("KMP_AFFINITY = scatter");

	int M	= A.size();
	if (B.size() != M)
	{
		Error_message_type	error;
		error.message		= "Problem in the size of Matrix for matrix_sub_vec_vec. PLEASE CHECK THE DIMENSION OF METRIX";
		error.module_name	= "MATRIX: OPPA_MATRIX.cpp";
		error.print_message();
	}

	vector<double>	C(M, 0.0);

#pragma omp parallel for num_threads(nt)
	for (int i = 0; i <= M-1; i++)	C[i]	= A[i] - B[i];

	return(C);
}



void block_tri_diagonal_decomposition(vector<vector<vector<double> > > &A,
									vector<vector<vector<double> > > &B,
									vector<vector<vector<double> > > &C, int nt)
{
	kmp_set_defaults("KMP_AFFINITY = scatter");

	int NB	= A.size();
	if (B.size() != NB || C.size() != NB)
	{
		Error_message_type	error;
		error.message		= "Problem in the size of Matrix for Tri-diagonal matrix decomposition. PLEASE CHECK THE DIMENSION OF METRIX";
		error.module_name	= "MATRIX: OPPA_MATRIX.cpp";
		error.print_message();
	}

	// Step 1: Error check
#pragma omp parallel for num_threads(nt)
	for (int i = 0; i <= NB-1; i++)
	{
		check_matrix_2D_square(A[i]);
		check_matrix_2D_square(B[i]);
		check_matrix_2D_square(C[i]);
	}



	// LU Decomposition for Block-tridiagonal matrix
	// L = |B1 0   0 ... 0|
	//   = |A2 B2  0 ... 0|
	//   = |   A3 B3 ... 0|
	//
	// U = |I C1  0  ... 0|
	//   = |0  I  C2 ... 0|
	//   = |0  0   I  C3 0|

	// Step 2: Get C_1
	vector< vector <double> > temp1 = matrix_inv(B[0]);
	vector< vector <double> > temp2	= matrix_multi_mat_mat(temp1, C[0], nt);;
	C[0]	= temp2;

	// Step 3: for other boundaries
	for (int i = 1; i <= NB-1; i++)
	{
		temp1	= matrix_multi_mat_mat(A[i], C[i-1], nt);
		temp2	= matrix_sub_mat_mat(B[i], temp1, nt);
		B[i]	= temp2;

		temp1	= matrix_inv(B[i]);
		temp2	= matrix_multi_mat_mat(temp1, C[i], nt);
		C[i]	= temp2;
	}
}


void block_tri_diagonal_solve(vector<vector<vector<double> > > &A,
								vector<vector<vector<double> > > &B,
								vector<vector<vector<double> > > &C,
								vector<vector <double> > &R,
								vector<vector <double> > &X,
								int nt)
{
	kmp_set_defaults("KMP_AFFINITY = scatter");

	int NB	= A.size();
	if (R.size() != NB)
	{
		Error_message_type	error;
		error.message		= "Problem in the size of Matrix for Tri-diagonal matrix solver. PLEASE CHECK THE DIMENSION OF METRIX";
		error.module_name	= "MATRIX: OPPA_MATRIX.cpp";
		error.print_message();
	}

	// Step 1:
	vector< vector <double> >	temp1 = matrix_inv(B[0]);
	vector< double> 			temp2 = matrix_multi_mat_vec(temp1, R[0]);
	X[0]	= temp2;

	// Step 2:
	vector< double> 			temp3;
	for (int i = 1; i <= NB-1; i++)
	{
		temp2 = matrix_multi_mat_vec(A[i], X[i-1]);
		temp3 = matrix_sub_vec_vec(R[i], temp2, nt);
		temp1 = matrix_inv(B[i]);

		temp2 = matrix_multi_mat_vec(temp1, temp3);
		X[i]	= temp2;
	}


	// Step 3:
	for (int i = NB-2; i >= 0; i--)
	{
		temp2 = matrix_multi_mat_vec(C[i], X[i+1]);
		temp3 = matrix_sub_vec_vec(X[i], temp3, nt);
		X[i] 	= temp3;;
	}
}


void block_tri_diagonal_solve_advance(vector<vector<vector<double> > > &A,
									vector<vector<vector<double> > > &B_inv,
									vector<vector<vector<double> > > &C,
									vector<vector <double> > &R,
									vector<vector <double> > &X,
									int nt)
{
	kmp_set_defaults("KMP_AFFINITY = scatter");

	int NB	= A.size();
	if (R.size() != NB)
	{
		Error_message_type	error;
		error.message		= "Problem in the size of Matrix for Tri-diagonal matrix solver. PLEASE CHECK THE DIMENSION OF METRIX";
		error.module_name	= "MATRIX: OPPA_MATRIX.cpp";
		error.print_message();
	}

	// Step 1:
	vector< double> 			temp2 = matrix_multi_mat_vec(B_inv[0], R[0]);
	X[0]	= temp2;

	// Step 2:
	vector< double> 			temp3;
	for (int i = 1; i <= NB-1; i++)
	{
		temp2 = matrix_multi_mat_vec(A[i], X[i-1]);
		temp3 = matrix_sub_vec_vec(R[i], temp2, nt);
		temp2 = matrix_multi_mat_vec(B_inv[i], temp3);

		X[i]	= temp2;
	}


	// Step 3:
	for (int i = NB-2; i >= 0; i--)
	{
		temp2 = matrix_multi_mat_vec(C[i], X[i+1]);
		temp3 = matrix_sub_vec_vec(X[i], temp2, nt);
		X[i] 	= temp3;;
	}
}



void Assemble_block_matrix(vector<vector<double> > &A,	vector<vector<double> > &M, int i, int j, int NVAR_I, int NVAR_J)
{
	int start_i	= NVAR_I * i;
	int start_j	= NVAR_J * j;


	for (int l = 0; l <= NVAR_I-1; l++)
		for (int m = 0; m <= NVAR_J-1; m++)
			M[start_i+l][start_j+m]	+= A[l][m];
}

