/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 10, 2014
 *      			Author: Minkwan Kim
 *
 * FMatrix.cpp
 * 			-  
 *  
 */

#include "../include/matrix.hpp"


void matrix_multi(vector < vector< double > > &A, vector < vector< double > > &B, vector < vector< double > > &C)
{
	int M	= A.size();

	for (int i = 0; i <= M-1; i++)
	{
		for (int j = 0; j <= M-1; j++)
		{
			double aux = 0.0;
			for (int k = 0; k <= M-1; k++)	aux += A[i][k]*B[k][j];

			C[i][j]	= aux;
		}
	}
}

// Constructors
matrix::matrix()
{
	N 			= 0;
	M 			= 0;
	elements	= NULL;
}

// Deconstructor
matrix::~matrix()
{
	if (M*N > 0)
	{
		DELETE_MAT<double>(elements, M, N);
		N	= 0;
		M	= 0;
	}

	elements = NULL;
}

matrix::matrix(int m, int n)
{
	M 			= m;
	N 			= n;
	elements	= NEW_MAT<double>(M, N);
}

matrix::matrix(const matrix &other)
{
	int i, j;
	N 			= other.N;
	M 			= other.M;
	elements	= NEW_MAT<double>(M, N);

	for (i = 0; i <= M-1; i++)
	{
		for (j = 0; j <= N-1; j++)
			elements[i][j] = other.elements[i][j];
	}
}

void matrix::destruct()
{
	if (M*N > 0)
	{
		DELETE_MAT<double>(elements, M, N);
		N	= 0;
		M	= 0;
	}
	elements = NULL;
}




//////////////////////////////////
/* Member function declarations */
//////////////////////////////////
//::M-01 - size()
//		 - ASSIGNING THE SIZE
void matrix::size(int m, int n)
{
	if ( m <= 0 || n <= 0)
	{
		throw Exception("[ERROR][matrix_MK]: Size of matrix should be unsigned integer!");
	}
	else
	{
		M			= m;
		N 			= n;
		elements	= NEW_MAT<double>(M, N);
	}
}


//::M-02 - diag()
// 		 - DIAGONAL MATRIX
void matrix::diag()
{
	int i, j;

	if (N*M == 0)
	{
		throw Exception("[ERROR][matrix_MK]-[diag()]: Size is not defined yet.");
	}
	else
	{
		if (M != N)
		{
			throw Exception("[ERROR][matrix_MK]-[diag()]: Diagonal natrix function is only applicable for a square matrix.");
		}
		else
		{
			for (i = 0; i <= M-1; i++)
			{
				for (j = 0; j <= N-1; j++)
				{
					if (i == j)	elements[i][j] = 1.0;
					else		elements[i][j] = 0.0;
				}
			}
		}
	}
}

void matrix::diag(int n)
{
	int i, j;

	if (N*M == 0)
	{
		size(n, n);
	}
	else
	{
		destruct();
		size(n, n);
	}

	diag();
}


//::M-03 - ones()
// 		 - ONES FUNCTIONS
void matrix::ones()
{
	int i, j;
	if (N*M == 0)
	{
		Exception("[ERROR][matrix_MK]: Size is not defined yet.");
	}
	else
	{
		for (i = 0; i <= M-1; i++)
			for (j = 0; j <= N-1; j++)
				elements[i][j] = 1.0;
	}
}

void matrix::ones(int n, int m)
{
	int i, j;

	if (N*M == 0)
	{
		size(n, n);
	}
	else
	{
		if (N != n || M != m)
		{
			destruct();
			size(n, m);
		}
	}

	ones();
}


//::M-04 - zeros()
// 		 - ZEROS FUNCTIONS
void matrix::zeros()
{
	int i, j;
	if (N*M == 0)
	{
		Exception("[ERROR][matrix_MK]: Size is not defined yet.");
	}
	else
	{
		for (i = 0; i <= M-1; i++)
			for (j = 0; j <= N-1; j++)
				elements[i][j] = 0.0;
	}
}

void matrix::zeros(int n, int m)
{
	int i, j;

	if (N*M == 0)
	{
		size(n, n);
	}
	else
	{
		if (N != n || M != m)
		{
			destruct();
			size(n, m);
		}
	}

	zeros();
}


//::M-05 - OPERATORS
double& matrix::operator()(const int i, const int j)
{
	if (i < 0 || i > M)
		throw Exception("[ERROR][matrix_MK]:  Subscript indices must either be real positive integers or logicals and it should be with a matrix size.");
	else if (j < 0 || j > N)
		throw Exception("[ERROR][matrix_MK]:  Subscript indices must either be real positive integers or logicals and it should be with a matrix size.");
	else
		return (elements[i][j]);
}

matrix& matrix::operator= (const matrix& a)
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


//::M-06 - ADD / SUBRACT / MULTIPLY
//		 - add
matrix& matrix::add(double v)
{
	int i, j;

	if (N*M == 0)
	{
		Exception("[ERROR][matrix_MK]: Size is not defined yet.");
	}
	else
	{
		for (i = 0; i <= M-1; i++)
			for (j = 0; j <= N-1; j++)
				elements[i][j] += v;
	}

	return *this;
}

matrix& matrix::sub(double v)
{
	return add(-v);
}

matrix& matrix::mul(double v)
{
	int i, j;

	if (N*M == 0)
	{
		Exception("[ERROR][matrix_MK]: Size is not defined yet.");
	}
	else
	{
		for (i = 0; i <= M-1; i++)
			for (j = 0; j <= N-1; j++)
				elements[i][j] *= v;
	}

	return *this;
}


//::M-07 - minor
//		 - MINOR FUNCTION
matrix matrix::minor_mk(const int r, const int c)
{
	int i, j;
	int res_r, res_c;

	matrix res;
	double temp;

	if (r >= 0 && r <= M-1 && c >= 0 && c <= N-1)
	{
		res.size(M-1, N-1);

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
						temp = elements[i][j];
						res.asg(res_r, res_c, temp);
						res_c++;
					}
				}
				res_r++;
			}
		}
	}

	return res;
}


//::M-08 - asg()
// 		 - ASSIGN DATA
void matrix::asg(int i, int j, double data)
{
	if (i < 0 || j < 0)
		Exception("[ERROR][matrix_MK]:  Subscript indices must either be real positive integers or logicals.");
	else
		elements[i][j] = data;
}



//////////////////////////////////
/* Global functions definitions */
//////////////////////////////////
// F-01: Operators
matrix operator +(const matrix& a, const matrix& b)
{
	int i, j;
	matrix res;

	if (a.M == b.M && a.N == b.N)
	{
		res.size(a.M, a.N);
		for (i = 0; i <= a.M-1; i++)
			for (j = 0;  j <= a.N-1; j++)
				res.elements[i][j] = a.elements[i][j] + b.elements[i][j];
	}
	else
	{
		Exception("[ERROR][matrix_MK]: Dimensions do not match.");
	}

	return res;
}

matrix operator+ (const matrix& A, const double B)
{
	matrix res = A;
	res.add(B);
	return res;
}

matrix operator+ (const double B, const matrix& A)
{
	matrix res = A;
	res.add(B);
	return res;
}

matrix operator- (const matrix& a, const matrix& b)
{
	int i, j;
	matrix res;

	if (a.M == b.M && a.N == b.N)
	{
		res.size(a.M, a.N);
		for (i = 0; i <= a.M-1; i++)
			for (j = 0;  j <= a.N-1; j++)
				res.elements[i][j] = a.elements[i][j] - b.elements[i][j];
	}
	else
	{
		Exception("[ERROR][matrix_MK]: Dimensions do not match.");
	}

	return res;
}

matrix operator- (const matrix& A, const double B)
{
	matrix res = A;
	res.sub(B);
	return res;
}

matrix operator- (const double B, const matrix& A)
{
	matrix res = A;
	res.sub(B);
	return res;
}

matrix operator- (const matrix& a)
{
	int i;
	matrix res;
	res.size(a.M, a.N);
	res = -1 * a;
	return res;
}

matrix operator* (const matrix& a, const matrix& b)
{
	int i, j, k;
	double temp;

	matrix res;
	if (a.N == b.M)
	{
		res.size(a.M, b.N);
		for (i = 0; i <= a.M-1; i++)
		{
			for (j = 0; j <= b.N-1; j++)
			{
				temp = 0;
				for (k = 0; k <= a.N-1; k++)
				{
					temp += a.elements[i][k]*b.elements[k][j];
				}
				res.elements[i][j] = temp;
			}
		}
		//res.size(a.M, b.N);
		//cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.M, a.N, b.N, 1.0, a.elements, a.N, a, ldb, 0.0, C, ldc);
	}
	else
	{
		Exception("[ERROR][matrix_MK]: Dimensions do not match.");
	}
	return res;
}


matrix operator* (const double B, const matrix& A)
{
	matrix res = A;
	res.mul(B);
	return res;
}


// F-02: ZEROS
matrix zeros(int i, int j)
{
	matrix A;
	A.zeros(i, j);
	return A;
}


// F-03: ONES
matrix ones(int i, int j)
{
	matrix A;
	A.ones(i, j);
	return A;
}


// F-04:: DIAGONAL MATRIX
matrix diag(int i)
{
	matrix A;
	A.zeros(i, i);
	for(int n = 1; n <= A.N; n++) A(n,n) = 1;

	return A;
}


// F-05:: Determinant
double det(matrix& A)
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
			d = A(0, 0);
		}
		else if (row == 2)
		{
			d = A(0,0) * A(1,1) - A(0,1) * A(1,0);
		}
		else
		{
			for (c = 0; c <= col-1; c++)
			{
				matrix M	= A.minor_mk(0, c);
				temp1 		= pow(-1.0, c);
				temp2 		= det(M);

				d += temp1 * A(0, c) * temp2;  // faster than with pow()
			}
		}
	}

	return d;
}


// F-06:: CONFACTOR
matrix confactor(matrix& a)
{
	int i, j;
	double d;
	matrix res;
	matrix ap;


	if (a.M == a.N && a.M != 0)
	{
		res.size(a.M, a.M);

		for (i = 0; i <= a.M-1; i++)
		{
			for (j = 0; j <= a.M-1; j++)
			{
				ap = a.minor_mk(i, j);
				d = det(ap);
				res(i,j) = pow(-1.0, i+j) * d;
			}
		}
	}
	else
	{
		Exception("[ERROR][matrix_MK]: Confactor matrix can be only calculated for square matrix");
	}

	ap.destruct();
	return res;
}


// F-07:: ADJOINT
matrix adjoint(matrix& a)
{
	int i, j;
	double temp;
	matrix res;
	matrix ap;

	ap = confactor(a);
	res.size(a.M, a.M);
	for (i = 0; i <= a.M-1; i++)
	{
		for (j = 0; j <= a.M-1; j++)
		{
			temp = ap(j,i);
			res(i,j) = temp;
		}
	}
	ap.destruct();
	return res;
}


// F-08:: INVERSE MATRIX
matrix inv(matrix& a)
{
	//int s;
	double d;
	matrix res;

	d = det(a);

	if (a.M == 1)		// 1 X 1 MATRIX
	{
		res.size(1,1);
		res(0,0) = 1/a(0,0);
	}
	else if (a.M == 2)	// 2 X 2 MATRIX
	{
		res.size(2,2);

		res(0,0) =  a(1,1);
		res(0,1) = -a(0,1);

		res(1,0) = -a(1,0);
		res(1,1) =  a(0,0);

		res = (1.0/d) * res;
	}
	else				// MORE THAN 3 X 3 MATRIX
	{
		// IT USES ADJOINT METHOD
		matrix ai;
		ai 	= adjoint(a);
		res = (1.0/d) * ai;
		ai.destruct();
		return res;
	}

	return res;
}







/*
 * ================================================
 * 		CMatrix class (using Intel MKL)
 * 											- Ver 1.0
 * ================================================
 */
// A. Constructors
CMatrix::CMatrix()
{
	M 			= 0;
	N 			= 0;
	elements	= NULL;
}

CMatrix::CMatrix(int m, int n)
{
	M 			= m;
	N 			= n;
	elements	= new double [m*n];
}

CMatrix::CMatrix(const CMatrix &other)
{
	int i;

	N 			= other.N;
	M 			= other.M;
	elements	= new double [M*N];

	for (i = 0; i <= M*N-1; i++)	elements[i] = other.elements[i];
}

CMatrix::CMatrix(const matrix &other)
{
	int i, j, k;
	N 			= other.N;
	M 			= other.M;
	elements	= new double [M*N];

	k = 0;
	for (i = 0; i <= M-1; i++)
	{
		for (j = 0; j <= N-1; j++)
		{
			elements[k] = other.elements[i][j];
			k++;
		}
	}
}


// B. Deconstructors
CMatrix::~CMatrix()
{
	if (M*N > 0)
	{
		delete [] elements;
		N	= 0;
		M	= 0;
	}

	elements = NULL;
}

void CMatrix::destruct()
{
	if (M*N > 0)
	{
		delete [] elements;
		N	= 0;
		M	= 0;
	}

	elements = NULL;
}



// Member functions
//
//::M-01 - size()
//		 - ASSIGNING THE SIZE
void CMatrix::size(int m, int n)
{
	if ( m <= 0 || n <= 0)
	{
		throw Exception("[ERROR][matrix_MK]: Size of matrix should be unsigned integer!");
	}
	else
	{
		M			= m;
		N 			= n;
		elements	= new double [M*N];
	}
}

//::M-02 - asg()
// 		 - ASSIGNING DATA
void CMatrix::asg(int i, int j, double data)
{
	int k	= Mat_index_2D_to_1D(i, j, M, N);
	elements[k] = data;
}


//::M-03 - val()
// 		 - Value of matrix
double CMatrix::val(int i, int j)
{
	int k	= Mat_index_2D_to_1D(i, j, M, N);
	return (elements[k]);
}


//::M-04 - diag()
// 		 - DIAGONAL MATRIX
void CMatrix::diag()
{
	int i, j, k;

	if (N*M == 0)	throw Exception("[ERROR][matrix_MK]-[diag()]: Size is not defined yet.");
	if (M != N)		throw Exception("[ERROR][matrix_MK]-[diag()]: Diagonal natrix function is only applicable for a square matrix.");


	for (i = 1; i <= M; i++)
	{
		for (j = 1; j <= N; j++)
		{
			k	= Mat_index_2D_to_1D(i, j, M, N);
			if (i == j)	elements[k] = 1.0;
			else		elements[k] = 0.0;
		}
	}
}

void CMatrix::diag(int n)
{

	if (N*M == 0)
	{
		size(n, n);
	}
	else
	{
		destruct();
		size(n, n);
	}

	diag();
}

//::M-05 - ones()
// 		 - ONES
void CMatrix::ones()
{
	int i, j, k;

	if (N*M == 0)
		throw Exception("[ERROR][matrix_MK]-[diag()]: Size is not defined yet.");

	for (i = 1; i <= M; i++)
		for (j = 1; j <= N; j++)
			asg(i, j, 1.0);
}

void CMatrix::ones(int m, int n)
{
	if (N*M == 0)
	{
		size(m, n);
	}
	else
	{
		if (N != n || M != m)
		{
			destruct();
			size(m, n);
		}
	}

	ones();
}

//::M-06 - zeros()
// 		 - ZEROS
void CMatrix::zeros()
{
	int i, j, k;

	if (N*M == 0)
		throw Exception("[ERROR][matrix_MK]-[diag()]: Size is not defined yet.");

	for (i = 1; i <= M; i++)
		for (j = 1; j <= N; j++)
			asg(i, j, 0.0);
}

void CMatrix::zeros(int m, int n)
{
	if (N*M == 0)
	{
		size(m, n);
	}
	else
	{
		if (N != n || M != m)
		{
			destruct();
			size(m, n);
		}
	}

	zeros();
}

//::M-07 - OPERATERS
double& CMatrix::operator()(const int i, const int j)
{
	if (i <= 0 || i > M)
	{
		throw Exception("[ERROR][matrix_MK]:  Subscript indices must either be real positive integers or logicals and it should be with a matrix size.");
	}
	else if (j <= 0 || j > N)
	{
		throw Exception("[ERROR][matrix_MK]:  Subscript indices must either be real positive integers or logicals and it should be with a matrix size.");
	}

	int k	= Mat_index_2D_to_1D(i, j, M, N);
	return (elements[k]);
}

CMatrix& CMatrix::operator= (const CMatrix& a)
{
	int i;
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
		for (i = 0; i <= N*M-1; i++)
			elements[i] = a.elements[i];
	}
	else
	{
		Exception("[ERROR][matrix_MK]: Two matrix should have same size.");
	}

	return *this;
}

CMatrix& CMatrix::operator= (const matrix& a)
{
	int i, j, k;
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
		for (i = 1; i <= M; i++)
		{
			for (j = 1; j <= N; j++)
			{
				asg(i, j, a.elements[i][j]);
			}
		}
	}
	else
	{
		Exception("[ERROR][matrix_MK]: Two matrix should have same size.");
	}

	return *this;
}


//::M-08 - ADD / SUBRACT / MULTIPLY
//		 - add
CMatrix& CMatrix::add(double v)
{
	int i;

	if (N*M == 0)
		Exception("[ERROR][matrix_MK]: Size is not defined yet.");

	for (i = 0; i <= M*N-1; i++)	elements[i] += v;

	return *this;
}

CMatrix& CMatrix::sub(double v)
{
	return add(-v);
}

CMatrix& CMatrix::mul(double v)
{
	int i;

	if (N*M == 0)
		Exception("[ERROR][matrix_MK]: Size is not defined yet.");

	for (i = 0; i <= M*N-1; i++)	elements[i] *= v;

	return *this;
}

//::M-07 - minor
//		 - MINOR FUNCTION
CMatrix CMatrix::minor_mk(const int r, const int c)
{
	int i, j;
	int res_r, res_c;

	CMatrix	res;

	if (r >= 1 && r <= M && c >= 1 && c <= N)
	{
		res.size(M-1, N-1);		// Create new matrix smaller size

		res_r = 0;
		for (i = 1; i <= M; i++)
		{
			if(i != r)
			{
				res_c = 0;
				for (j = 1; j <= N; j++)
				{
					if(j != c)
					{
						res.asg(res_r, res_c, val(i, j));
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
 * CVector class (using Intel MKL)
 * 		- Ver 0.0
 */
// A. Constructors
CVector::CVector()
{
	M 			= 0;
	elements	= NULL;
}

CVector::CVector(int m)
{
	M 			= m;
	elements	= new double [m];
}

CVector::CVector(const CVector &other)
{
	int i;

	M 			= other.M;
	elements	= new double [M];

	for (i = 0; i <= M-1; i++)	elements[i] = other.elements[i];
}

CVector::CVector(int n, double *a)
{
	int i;
	M 			= n;
	elements	= new double [M];

	for (i = 0; i <= M-1; i++)
	{
		elements[i] = a[i];
	}
}

// B. Deconstructors
CVector::~CVector()
{
	if (M > 0)
	{
		delete [] elements;
		M	= 0;
	}

	elements = NULL;
}

void CVector::destruct()
{
	if (M > 0)
	{
		delete [] elements;
		M	= 0;
	}

	elements = NULL;
}
















//////////////////////////////////
/* Global functions definitions */
//////////////////////////////////
// CF-01: Operators
CMatrix operator +(const CMatrix& A, const CMatrix& B)
{
	int i, j;
	CMatrix res;

	if (A.M == A.M && A.N == B.N)
	{
		res.size(A.M, B.N);
		for (i = 0; i <= A.M*A.N-1; i++)
			res.elements[i] = A.elements[i] + B.elements[i];
	}
	else
	{
		Exception("[ERROR][matrix_MK]: Dimensions do not match.");
	}

	return res;
}

CMatrix operator+ (const CMatrix& A, const double B)
{
	CMatrix res = A;
	res.add(B);
	return res;
}

CMatrix operator+ (const double B, const CMatrix& A)
{
	CMatrix res = A;
	res.add(B);
	return res;
}


CMatrix operator- (const CMatrix& A, const CMatrix& B)
{
	int i, j;
	CMatrix res;

	if (A.M == B.M && A.N == B.N)
	{
		res.size(A.M, B.N);
		for (i = 0; i <= A.M*A.N-1; i++)
			res.elements[i] = A.elements[i] - B.elements[i];
	}
	else
	{
		Exception("[ERROR][matrix_MK]: Dimensions do not match.");
	}

	return res;
}

CMatrix operator- (const CMatrix& A, const double B)
{
	CMatrix res = A;
	res.sub(B);
	return res;
}

CMatrix operator- (const double B, const CMatrix& A)
{
	CMatrix res = A;
	res.sub(B);
	return res;
}

CMatrix operator- (const CMatrix& A)
{
	CMatrix res;
	res.size(A.M, A.N);
	res = -1.0 * A;
	return res;
}

CMatrix operator* (const CMatrix& A, const CMatrix& B)
{
	int m, n, k;
	double alpha, beta;

	CMatrix res;
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


CMatrix operator* (const double B, const CMatrix& A)
{
	CMatrix res = A;
	res.mul(B);
	return res;
}


// CF-02: ZEROS
CMatrix ZEROS(int i, int j)
{
	CMatrix A;
	A.zeros(i, j);
	return A;
}

CMatrix ONES(int i, int j)
{
	CMatrix A;
	A.ones(i, j);
	return A;
}

CMatrix DIAG(int i)
{
	CMatrix A;
	A.diag(i);

	return A;
}

// CF-03:: INVERSE MATRIX
CMatrix INV(CMatrix& a)
{
	int 	i, j;
	int 	info;
	int 	*ipiv;
	int		m;
	CMatrix A_inv;

	if (a.M != a.N)
		Exception("[ERROR][matrix_MK]: Square matrix is required to calculate inverse matrix");

	m	= a.M;
	A_inv.size(m, m);

	info			= 0;
	ipiv			= new int [a.M+1];
	int 	lwork	= A_inv.M*A_inv.M;
	double *work	= new double [lwork];
	double *element = new double [m*m];

	// Store element in Column-major order
	for (i = 0; i <= m-1; i++)
		for (j = 0; j <= m-1; j++)
			element[i + j*m] = a(i+1,j+1);

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
double DET(CMatrix& A)
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
			d = A.val(1, 1);
		}
		else if (row == 2)
		{
			d = A.val(1,1) * A.val(2,2) - A.val(1,2) * A.val(2,1);
		}
		else
		{
			for (c = 1; c <= col; c++)
			{
				CMatrix M	= A.minor_mk(1, c);
				temp1 		= pow(-1.0, c);
				temp2 		= DET(M);

				d += temp1 * A.val(1, c) * temp2;  // faster than with pow()
			}
		}
	}

	return d;
}


// CF-05:: Block Tridiagonal matric decomposition
void block_tri_diagonal_decomp(CMatrix *A, CMatrix *B, CMatrix *C, int NB)
{
	int 	i;
	int 	N;
	int		error_index;
	bool	error_flag;

	// Step 1: Error check
	error_flag = false;
	for (i = 1; i <= NB; i++)
	{
		if (A[i].M != A[i].N)
		{
			error_flag	= true;
			error_index	= 1;
			break;
		}

		if (B[i].M != B[i].N)
		{
			error_flag	= true;
			error_index	= 2;
			break;
		}

		if (C[i].M != C[i].N)
		{
			error_flag	= true;
			error_index	= 3;
			break;
		}
	}

	if (error_flag == true)
	{
		Error_message_type	error;
		error.location_primary_name	="matrix.cpp";
		error.message = "Error in Block Tridiagonal Matrix solver: The Block matrices should be square!!";
		error.print_message();
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
	C[1]	= INV(B[1]) * C[1];

	// Step 3: for other boundaries
	for (i = 2; i <= NB; i++)
	{
		B[i]	= B[i] - A[i]*C[i-1];
		C[i]	= INV(B[i]) * C[i];
	}
}

void block_tri_diagonal_decomp2(vector<CMatrix> &A, vector<CMatrix> &B, vector<CMatrix> &C, int NB)
{
	int 	i;
	int 	N;
	int		error_index;
	bool	error_flag;

	// Step 1: Error check
	error_flag = false;
	for (i = 1; i <= NB; i++)
	{
		if (A[i].M != A[i].N)
		{
			error_flag	= true;
			error_index	= 1;
			break;
		}

		if (B[i].M != B[i].N)
		{
			error_flag	= true;
			error_index	= 2;
			break;
		}

		if (C[i].M != C[i].N)
		{
			error_flag	= true;
			error_index	= 3;
			break;
		}
	}

	if (error_flag == true)
	{
		Error_message_type	error;
		error.location_primary_name	="matrix.cpp";
		error.message = "Error in Block Tridiagonal Matrix solver: The Block matrices should be square!!";
		error.print_message();
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
	C[1]	= INV(B[1]) * C[1];

	// Step 3: for other boundaries
	for (i = 2; i <= NB; i++)
	{
		B[i]	= B[i] - A[i]*C[i-1];
		C[i]	= INV(B[i]) * C[i];
	}
}



void block_tri_diagonal_decomp3(vector<CMatrix> &A, vector<CMatrix> &B, vector<CMatrix> &C, int NB)
{
	int 	i;
	int 	N;
	int		error_index;
	bool	error_flag;

	// Step 1: Error check
	error_flag = false;
	for (i = 0; i <= NB-1; i++)
	{
		if (A[i].M != A[i].N)
		{
			error_flag	= true;
			error_index	= 1;
			break;
		}

		if (B[i].M != B[i].N)
		{
			error_flag	= true;
			error_index	= 2;
			break;
		}

		if (C[i].M != C[i].N)
		{
			error_flag	= true;
			error_index	= 3;
			break;
		}
	}

	if (error_flag == true)
	{
		Error_message_type	error;
		error.location_primary_name	="matrix.cpp";
		error.message = "Error in Block Tridiagonal Matrix solver: The Block matrices should be square!!";
		error.print_message();
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
	C[0]	= INV(B[0]) * C[0];

	// Step 3: for other boundaries
	for (i = 1; i <= NB-1; i++)
	{
		B[i]	= B[i] - A[i]*C[i-1];
		C[i]	= INV(B[i]) * C[i];
	}
}



// CF-06:: Block Tridiagonal matrices Solver
void block_tri_diagonal_solver(CMatrix *A, CMatrix *B, CMatrix *C, CMatrix *R, CMatrix *X, int NB)
{
	int 	i;
	int 	N;
	int		error_index;
	bool	error_flag;


	// Step 1: Error check
	error_flag = false;
	for (i = 1; i <= NB; i++)
	{
		if (A[i].M != A[i].N)
		{
			error_flag	= true;
			error_index	= 1;
			break;
		}

		if (B[i].M != B[i].N)
		{
			error_flag	= true;
			error_index	= 2;
			break;
		}

		if (C[i].M != C[i].N)
		{
			error_flag	= true;
			error_index	= 3;
			break;
		}
	}

	if (error_flag == true)
	{
		Error_message_type	error;
		error.location_primary_name	="matrix.cpp";
		error.message = "Error in Block Tridiagonal Matrix solver: The Block matrices should be square!!";
		error.print_message();
	}

	N	= A[1].M;
	for (i = 1; i <= NB; i++)
	{
		if (R[i].N != 1)
		{
			error_flag	= true;
			error_index	= 1;
			break;
		}

		if (R[i].M != N)
		{
			error_flag	= true;
			error_index	= 2;
			break;
		}

		if (X[i].N != 1)
		{
			error_flag	= true;
			error_index	= 3;
			break;
		}

		if (X[i].M != N)
		{
			error_flag	= true;
			error_index	= 4;
			break;
		}
	}

	if (error_flag == true)
	{
		cout << "Error in Block Tridiagonal Matrix solver:";
		switch(error_index)
		{
		case 1:
			cout <<" R should be vector (not matrix)" << endl;
			break;
		case 2:
			cout <<" The size of vector R should be same as A" << endl;
			break;
		case 3:
			cout <<" X should be vector (not matrix)" << endl;
			break;
		case 4:
			cout <<" The size of vector X should be same as A" << endl;
			break;
		}

		Error_message_type	error;
		error.location_primary_name	="matrix.cpp";
		error.message = "Please check block matrices!";
		error.print_message();
	}


	// Step 1:
	X[1]	= B[1] * R[1];

	// Step 2:
	for (i = 2; i <= NB; i++)
		X[i]	= B[i] * (R[i] - A[i]*X[i-1]);

	// Step 3:
	for (i = NB-1; i >= 1; i--)
	{
		X[i] 	= X[i] - C[i]*X[i+1];
	}
}

void block_tri_diagonal_solver2(vector<CMatrix> &A, vector<CMatrix> &B, vector<CMatrix> &C, vector<CMatrix> &R, vector<CMatrix> &X, int NB)
{
	int 	i;
	int 	N;
	int		error_index;
	bool	error_flag;


	// Step 1: Error check
	error_flag = false;
	for (i = 1; i <= NB; i++)
	{
		if (A[i].M != A[i].N)
		{
			error_flag	= true;
			error_index	= 1;
			break;
		}

		if (B[i].M != B[i].N)
		{
			error_flag	= true;
			error_index	= 2;
			break;
		}

		if (C[i].M != C[i].N)
		{
			error_flag	= true;
			error_index	= 3;
			break;
		}
	}

	if (error_flag == true)
	{
		Error_message_type	error;
		error.location_primary_name	="matrix.cpp";
		error.message = "Error in Block Tridiagonal Matrix solver: The Block matrices should be square!!";
		error.print_message();
	}

	N	= A[1].M;
	for (i = 1; i <= NB; i++)
	{
		if (R[i].N != 1)
		{
			error_flag	= true;
			error_index	= 1;
			break;
		}

		if (R[i].M != N)
		{
			error_flag	= true;
			error_index	= 2;
			break;
		}

		if (X[i].N != 1)
		{
			error_flag	= true;
			error_index	= 3;
			break;
		}

		if (X[i].M != N)
		{
			error_flag	= true;
			error_index	= 4;
			break;
		}
	}

	if (error_flag == true)
	{
		cout << "Error in Block Tridiagonal Matrix solver:";
		switch(error_index)
		{
		case 1:
			cout <<" R should be vector (not matrix)" << endl;
			break;
		case 2:
			cout <<" The size of vector R should be same as A" << endl;
			break;
		case 3:
			cout <<" X should be vector (not matrix)" << endl;
			break;
		case 4:
			cout <<" The size of vector X should be same as A" << endl;
			break;
		}

		Error_message_type	error;
		error.location_primary_name	="matrix.cpp";
		error.message = "Please check block matrices!!";
		error.print_message();
	}


	// Step 1:
	X[1]	= B[1] * R[1];

	// Step 2:
	for (i = 2; i <= NB; i++)
		X[i]	= B[i] * (R[i] - A[i]*X[i-1]);

	// Step 3:
	for (i = NB-1; i >= 1; i--)
	{
		X[i] 	= X[i] - C[i]*X[i+1];
	}
}


void block_tri_diagonal_solver3(vector<CMatrix> &A, vector<CMatrix> &B, vector<CMatrix> &C, vector<CMatrix> &R, vector<CMatrix> &X, int NB)
{
	int 	i;
	int 	N;
	int		error_index;
	bool	error_flag;


	// Step 1: Error check
	error_flag = false;
	for (i = 0; i <= NB-1; i++)
	{
		if (A[i].M != A[i].N)
		{
			error_flag	= true;
			error_index	= 1;
			break;
		}

		if (B[i].M != B[i].N)
		{
			error_flag	= true;
			error_index	= 2;
			break;
		}

		if (C[i].M != C[i].N)
		{
			error_flag	= true;
			error_index	= 3;
			break;
		}
	}

	if (error_flag == true)
	{
		Error_message_type	error;
		error.location_primary_name	="matrix.cpp";
		error.message = "Error in Block Tridiagonal Matrix solver: The Block matrices should be square!!";
		error.print_message();
	}

	N	= A[0].M;
	for (i = 0; i <= NB-1; i++)
	{
		if (R[i].N != 1)
		{
			error_flag	= true;
			error_index	= 1;
			break;
		}

		if (R[i].M != N)
		{
			error_flag	= true;
			error_index	= 2;
			break;
		}

		if (X[i].N != 1)
		{
			error_flag	= true;
			error_index	= 3;
			break;
		}

		if (X[i].M != N)
		{
			error_flag	= true;
			error_index	= 4;
			break;
		}
	}

	if (error_flag == true)
	{
		cout << "Error in Block Tridiagonal Matrix solver:";
		switch(error_index)
		{
		case 1:
			cout <<" R should be vector (not matrix)" << endl;
			break;
		case 2:
			cout <<" The size of vector R should be same as A" << endl;
			break;
		case 3:
			cout <<" X should be vector (not matrix)" << endl;
			break;
		case 4:
			cout <<" The size of vector X should be same as A" << endl;
			break;
		}

		Error_message_type	error;
		error.location_primary_name	="matrix.cpp";
		error.message = "Please check block matrices!";
		error.print_message();
	}


	// Step 1:
	X[0]	= B[0] * R[0];

	// Step 2:
	for (i = 1; i <= NB-1; i++)	X[i]	= B[i] * (R[i] - A[i]*X[i-1]);

	// Step 3:
	for (i = NB-2; i >= 0; i--)	X[i] 	= X[i] - C[i]*X[i+1];
}







// CF-07:: Delete CMatrix pointer
void DELETE_CMatrix_p(CMatrix *A, int N)
{
	int i;
	for (i = 0; i <= N-1; i++)
	{
		A[i].destruct();
	}

	delete [] A;
}



// CF-08:: Block Tridiagonal matrices Solver
void F_Block_tri_diagonal_solver(CMatrix *A, CMatrix *B, CMatrix *C, CMatrix *R, CMatrix *X, int NB)
{
	int 	i;
	int 	N;
	int		error_index;
	bool	error_flag;


	// Step 1: Error check
	error_flag	= false;
	N			= A[1].M;
	for (i = 1; i <= NB; i++)
	{
		if (R[i].N != 1)
		{
			error_flag	= true;
			error_index	= 1;
			break;
		}

		if (R[i].M != N)
		{
			error_flag	= true;
			error_index	= 2;
			break;
		}

		if (X[i].N != 1)
		{
			error_flag	= true;
			error_index	= 3;
			break;
		}

		if (X[i].M != N)
		{
			error_flag	= true;
			error_index	= 4;
			break;
		}
	}

	if (error_flag == true)
	{
		cout << "Error in Block Tridiagonal Matrix solver:";
		switch(error_index)
		{
		case 1:
			cout <<" R should be vector (not matrix)" << endl;
			break;
		case 2:
			cout <<" The size of vector R should be same as A" << endl;
			break;
		case 3:
			cout <<" X should be vector (not matrix)" << endl;
			break;
		case 4:
			cout <<" The size of vector X should be same as A" << endl;
			break;
		}

		Error_message_type	error;
		error.location_primary_name	="matrix.cpp";
		error.message = "Please check block matrices!";
		error.print_message();
	}


	// Step 2:: LU Factorization
	block_tri_diagonal_decomp(A, B, C, NB);

	// Step 3:: Forward substitution
	CMatrix *Y;
	Y	= new CMatrix [NB+1];
	for (i = 0; i <= NB; i++)	Y[i].zeros(N, 1);

	CMatrix aux;
	aux = INV(B[1])*B[1];

	Y[1]	= INV(B[1]) * R[1];
	for (i = 2;  i <= NB; i++)
		Y[i]	= INV(B[i]) * (R[i] - A[i]*Y[i-1]);


	// Step 4:: Backward substitution
	X[NB]	= Y[NB];
	for (i = NB-1; i >= 1; i--)
		X[i]= Y[i] - C[i] * X[i+1];


	// Step Final:: Delete allocated matrix
	DELETE_CMatrix_p(Y, NB+1);
}


// CF-08:: Block Tridiagonal matrices Solver
void F_Block_tri_diagonal_solver2(vector<CMatrix> &A, vector<CMatrix> B, vector<CMatrix> &C, vector<CMatrix> &R, vector<CMatrix> &X, int NB)
{
	int 	i;
	int 	N;
	int		error_index;
	bool	error_flag;


	// Step 1: Error check
	error_flag	= false;
	N			= A[1].M;
	for (i = 1; i <= NB; i++)
	{
		if (R[i].N != 1)
		{
			error_flag	= true;
			error_index	= 1;
			break;
		}

		if (R[i].M != N)
		{
			error_flag	= true;
			error_index	= 2;
			break;
		}

		if (X[i].N != 1)
		{
			error_flag	= true;
			error_index	= 3;
			break;
		}

		if (X[i].M != N)
		{
			error_flag	= true;
			error_index	= 4;
			break;
		}
	}

	if (error_flag == true)
	{
		cout << "Error in Block Tridiagonal Matrix solver:";
		switch(error_index)
		{
		case 1:
			cout <<" R should be vector (not matrix)" << endl;
			break;
		case 2:
			cout <<" The size of vector R should be same as A" << endl;
			break;
		case 3:
			cout <<" X should be vector (not matrix)" << endl;
			break;
		case 4:
			cout <<" The size of vector X should be same as A" << endl;
			break;
		}
		Error_message_type	error;
		error.location_primary_name	="matrix.cpp";
		error.message = "Please check block matrices!";
		error.print_message();
	}


	// Step 2:: LU Factorization
	block_tri_diagonal_decomp2(A, B, C, NB);

	// Step 3:: Forward substitution
	vector<CMatrix> Y(NB+1);
	for (i = 0; i <= NB; i++)	Y[i].zeros(N, 1);

	CMatrix aux;
	aux = INV(B[1])*B[1];

	Y[1]	= INV(B[1]) * R[1];
	for (i = 2;  i <= NB; i++)
		Y[i]	= INV(B[i]) * (R[i] - A[i]*Y[i-1]);


	// Step 4:: Backward substitution
	X[NB]	= Y[NB];
	for (i = NB-1; i >= 1; i--)
		X[i]= Y[i] - C[i] * X[i+1];

}


















vector< vector<double> > INV(vector< vector<double> >& a)
{
	int 	i, j;
	int 	info;
	int 	*ipiv;
	int 	m	= a.size();
	vector <vector <double> >	a_inv	= vector_2D(m, m, 0.0);

	info			= 0;
	ipiv			= new int [m+1];

	int 	lwork	= m*m;
	double *work	= new double [lwork];
	double *element = new double [m*m];

	// Store element in Column-major order
	for (i = 0; i <= m-1; i++)	for (j = 0; j <= m-1; j++)	element[i + j*m] = a[i][j];

	dgetrf(&m, &m, element, &m, ipiv, &info);
	dgetri(&m, element, &m, ipiv, work, &lwork, &info);

	int k;
	k	= 0;
	for (j = 0; j <= m-1; j++)
	{
		for (i = 0; i <= m-1; i++)
		{
			a_inv[i][j] = element[k];
			k++;
		}
	}

	delete[] ipiv;
	delete[] work;
	delete[] element;

	return (a_inv);
}




void Assemble_block_matrix2(CMatrix &A,	CMatrix &M, int i, int j, int NVAR_I, int NVAR_J)
{
	int start_i	= NVAR_I*(i-1);
	int start_j	= NVAR_J*(j-1);


	for (int l = 1; l <= NVAR_I; l++)
		for (int m = 1; m <= NVAR_J; m++)
			M(start_i+l, start_j+m)	+= A(l, m);
}


