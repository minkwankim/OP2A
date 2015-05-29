/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Mar 13, 2015
 *      			Author: Minkwan Kim
 *
 * OPPA_Matrix.hpp
 * 			-  
 *  
 */

#ifndef OPPA_MATRIX2_HPP_
#define OPPA_MATRIX2_HPP_

#include "mkl.h"
#include <cmath>
#include <sstream>
#include <iostream>
#include <vector>
#include "../../UTIL/include/OP2A_utilities.hpp"
#include "math_misc.hpp"


using namespace std;


void check_matrix_2D_square(vector<vector<double> > &A);
vector<vector<double> > matrix_minor(vector<vector<double> > &A, int r, int c);
double matrix_det(vector<vector<double> > &A);
vector<vector<double> > matrix_confactor(vector <vector<double> > &A);
vector<vector<double> > matrix_adjoint(vector<vector<double> >& A);
vector<vector<double> > matrix_inv(vector<vector<double> > &A);
vector<vector<double> > matrix_inv2(vector<vector<double> > &A);

vector<double> matrix_multi_mat_vec(vector<vector<double> > &A, vector<double> &X);
vector< vector<double> > matrix_multi_mat_mat(vector<vector<double> > &A, vector<vector<double> > &B, int nt);
vector< vector<double> > matrix_sub_mat_mat(vector<vector<double> > &A, vector<vector<double> > &B, int nt);
vector<double> matrix_sub_vec_vec(vector<double> &A, vector<double> &B, int nt);

void block_tri_diagonal_decomposition(vector<vector<vector<double> > > &A, vector<vector<vector<double> > > &B, vector<vector<vector<double> > > &C, int nt);
void block_tri_diagonal_solve(vector<vector<vector<double> > > &A, vector<vector<vector<double> > > &B, vector<vector<vector<double> > > &C, vector<vector <double> > &R, vector<vector <double> > &X, int nt);
void block_tri_diagonal_solve_advance(vector<vector<vector<double> > > &A, vector<vector<vector<double> > > &B_inv, vector<vector<vector<double> > > &C, vector<vector <double> > &R, vector<vector <double> > &X, int nt);
void Assemble_block_matrix(vector<vector<double> > &A,	vector<vector<double> > &M, int i, int j, int NVAR_I, int NVAR_J);


/* MATRIX CLASS */
class Matrix
{
public:
	vector < vector<double> > elements;
	int M;
	int N;

	// CONSTRACTOR & DECONSTRUCTOR
	Matrix();
	Matrix(int m, int n);
	Matrix(const Matrix &other);
	Matrix(vector< vector <double> > &other);



	~Matrix();
	void destruct();

	//////////////////////////////////
	// Member function declarations	//
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
	double&	operator()(const int i, const int j);
	Matrix&	operator= (const Matrix& a);
	Matrix&	operator= (vector< vector<double> >& a);

	//::M-06 - ADD / SUBRACT / MULTIPLY
	//		 - add
	Matrix& add(double v);
	Matrix& sub(double v);
	Matrix& mul(double v);

	//::M-07 - minor
	//		 - MINOR FUNCTION
	Matrix minor_mk(const int r, const int c);

};

Matrix	operator* (const Matrix& A, 	const Matrix& B);


#endif /* OPPA_MATRIX_HPP_ */
