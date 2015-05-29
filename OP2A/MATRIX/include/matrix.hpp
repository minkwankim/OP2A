/************************************************************************
			Open-source Multi-Physics Solver - ver. 0.0		
				
Copyright(c) 2013 MINKWAN KIM												
							
Initial Developed Date: Oct 29, 2013
                    by: Minkwan Kim

Last Modified Date: Oct 28, 2013
					by: Minkwan Kim


matrix.h
	- Library for MATRIX calculations
              
					Copyright 2013 MINKWAN KIM
*************************************************************************/

#ifndef _MATRIX_H
#define _MATRIX_H

#include "mkl.h"
#include <cmath>
#include <sstream>
#include <iostream>
#include <vector>

#include "../../UTIL/include/OP2A_utilities.hpp"
#include "CMatrix_basic.hpp"

using namespace std;



/* MATRIX CLASS */
class matrix
{
public:
	double** elements;
	int M;
	int N;
	
	// CONSTRACTOR & DECONSTRUCTOR
	matrix();
	~matrix();

	matrix(int n, int m);
	matrix(matrix const &other);
	void destruct();

	//////////////////////////////////
	/* Member function declarations */
	//////////////////////////////////
	//::M-01 - size()
	//		 - ASSIGNING THE SIZE
	void size(int m, int n);

	//::M-02 - diag()
	// 		 - DIAGONAL MATRIX
	void diag();
	void diag(int n);

	//::M-03 - ones()
	// 		 - ONES FUNCTIONS
	void ones();
	void ones(int i, int j);
	
	//::M-04 - zeros()
	// 		 - ZEROSS FUNCTIONS
	void zeros();
	void zeros(int i, int j);

	//::M-05 - OPERATERS
	double& operator()(const int i, const int j);
	matrix&	operator= (const matrix& a);

	//::M-06 - ADD / SUBRACT / MULTIPLY
	//		 - add
	matrix& add(double v);
	matrix& sub(double v);
	matrix& mul(double v);
	
	//::M-07 - minor
	//		 - MINOR FUNCTION
	matrix minor_mk(const int r, const int c);

	//::M-08 - asg()
	// 		 - ASSIGN DATA
	void asg(int i, int j, double data);
};


///////////////////////////////////
/* Global functions declarations */
///////////////////////////////////
// F-01: Operators
matrix	operator+ (const matrix& a, const matrix& b);
matrix	operator+ (const matrix& A, const double B);
matrix	operator+ (const double B, 	const matrix& A);
matrix	operator- (const matrix& a, const matrix& b);
matrix	operator- (const matrix& A, const double B);
matrix 	operator- (const double B, const matrix& A);
matrix	operator- (const matrix& a);
matrix	operator* (const matrix& a, const matrix& b);
matrix	operator* (const double B, const matrix& A);

// F-02: ZEROS
matrix zeros(int i, int j);

// F-03: ONES
matrix ones(int i, int j);

// F-04:: DIAGONAL MATRIX
matrix diag(int i);

// F-05:: Determinant
double det(matrix& A);

// F-06:: CONFACTOR
matrix confactor(matrix& a);

// F-07:: ADJOINT
matrix adjoint(matrix& a);

// F-08:: INVERSE MATRIX
matrix inv(matrix& a);







/* MATRIX CLASS */
class CMatrix
{
public:
	double* elements;
	int M;
	int N;

	// CONSTRACTOR & DECONSTRUCTOR
	CMatrix();
	CMatrix(int m, int n);
	CMatrix(const CMatrix &other);
	CMatrix(const matrix &other);

	~CMatrix();
	void destruct();



	//////////////////////////////////
	// Member function declarations	//
	//////////////////////////////////
	//::M-01 - size()
	//		 - ASSIGNING THE SIZE
	void size(int m, int n);

	//::M-02 - asg()
	// 		 - ASSIGN DATA
	void asg(int i, int j, double data);

	//::M-03 - val()
	// 		 - Value of matrix
	double val(int i, int j);

	//::M-04 - diag()
	// 		 - DIAGONAL MATRIX
	void diag();
	void diag(int n);

	//::M-05 - ones()
	// 		 - ONES FUNCTIONS
	void ones();
	void ones(int i, int j);

	//::M-06 - zeros()
	// 		 - ZEROSS FUNCTIONS
	void zeros();
	void zeros(int i, int j);

	//::M-07 - OPERATERS
	double& 	operator()(const int i, const int j);
	CMatrix&	operator= (const CMatrix& a);
	CMatrix&	operator= (const matrix& a);

	//::M-08 - ADD / SUBRACT / MULTIPLY
	//		 - add
	CMatrix& add(double v);
	CMatrix& sub(double v);
	CMatrix& mul(double v);

	//::M-08 - minor
	//		 - MINOR FUNCTION
	CMatrix minor_mk(const int r, const int c);

};



/* Vector CLASS */
class CVector
{
public:
	double* elements;
	int M;

	// CONSTRACTOR & DECONSTRUCTOR
	CVector();
	CVector(int m);
	CVector(const CVector &other);
	CVector(int n, double *a);

	~CVector();
	void destruct();

};










///////////////////////////////////
/* Global functions declarations */
///////////////////////////////////
// CF-01: Operators
CMatrix	operator+ (const CMatrix& A, 	const CMatrix& B);
CMatrix	operator+ (const CMatrix& A, 	const double B);
CMatrix	operator+ (const double B, 		const CMatrix& A);

CMatrix	operator- (const CMatrix& A, 	const CMatrix& B);
CMatrix	operator- (const CMatrix& A, 	const double B);
CMatrix operator- (const double B, 		const CMatrix& A);
CMatrix	operator- (const CMatrix& A);

CMatrix	operator* (const CMatrix& A, 	const CMatrix& B);
CMatrix	operator* (const double B,		const CMatrix& A);

// CF-02: ZEROS/ONES/DIAGONAL
CMatrix ZEROS(int i, int j);
CMatrix ONES(int i, int j);
CMatrix DIAG(int i);

// CF-03:: INVERSE MATRIX
CMatrix INV(CMatrix& a);
vector< vector<double> > INV(vector< vector<double> >& a);

// CF-04:: Determinant
double DET(CMatrix& A);

// CF-05:: Block Tridiagonal matrices decomposition
void block_tri_diagonal_decomp(CMatrix *A, CMatrix *B, CMatrix *C, int NB);
void block_tri_diagonal_decomp2(vector<CMatrix> &A, vector<CMatrix> &B, vector<CMatrix> &C, int NB);
void block_tri_diagonal_decomp3(vector<CMatrix> &A, vector<CMatrix> &B, vector<CMatrix> &C, int NB);

// CF-06:: Block Tridiagonal matrices Solver
void block_tri_diagonal_solver(CMatrix *A, CMatrix *B, CMatrix *C, CMatrix *R, CMatrix *X, int NB);
void block_tri_diagonal_solver2(vector<CMatrix> &A, vector<CMatrix> &B, vector<CMatrix> &C, vector<CMatrix> &R, vector<CMatrix> &X, int NB);
void block_tri_diagonal_solver3(vector<CMatrix> &A, vector<CMatrix> &B, vector<CMatrix> &C, vector<CMatrix> &R, vector<CMatrix> &X, int NB);

// CF-07:: Delete CMatrix pointer
void DELETE_CMatrix_p(CMatrix *A, int N);

// CF-08:: Block Tridiagonal matrices Solver
void F_Block_tri_diagonal_solver(CMatrix *A, CMatrix *B, CMatrix *C, CMatrix *R, CMatrix *X, int NB);
void F_Block_tri_diagonal_solver2(vector<CMatrix> &A, vector<CMatrix> B, vector<CMatrix> &C, vector<CMatrix> &R, vector<CMatrix> &X, int NB);


void Assemble_block_matrix2(CMatrix &A,	CMatrix &M, int i, int j, int NVAR_I, int NVAR_J);





#endif
