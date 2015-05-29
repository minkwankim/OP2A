/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Mar 23, 2015
 *      			Author: Minkwan Kim
 *
 * CMatrix.hpp
 * 			-  
 *  
 */

#ifndef CMATRIX2_HPP_
#define CMATRIX2_HPP_

#include "mkl.h"
#include <cmath>
#include <sstream>
#include <iostream>
#include <vector>


#include "../../UTIL/include/OP2A_utilities.hpp"
#include "./CMatrix_basic.hpp"


class CMatrix_ver2
{
public:
	double* elements;
	int M;
	int N;

	// CONSTRACTOR & DECONSTRUCTOR
	CMatrix_ver2();
	CMatrix_ver2(int m, int n);
	CMatrix_ver2(vector<vector<double> > &other);

	~CMatrix_ver2();



	//////////////////////////////////
	// Member function declarations	//
	//////////////////////////////////
	//::M-01 - size()
	//		 - ASSIGNING THE SIZE
	void size(int m, int n);

	//::M-01 - OPERATERS
	double&		operator() (const int i, const int j);
	CMatrix_ver2&	operator= (const CMatrix_ver2 &A);
	CMatrix_ver2&	operator= (const vector< vector<double> > &A);
	CMatrix_ver2&	operator+= (const double a);
	CMatrix_ver2&	operator+= (const CMatrix_ver2& A);
	CMatrix_ver2&	operator-= (const double a);
	CMatrix_ver2&	operator-= (const CMatrix_ver2& A);

	//::M-02 - Basic creation functions
	// 		 - DIAGONAL MATRIX
	void diag();
	void diag(int n);

	// 		 - ONES FUNCTIONS
	void ones();
	void ones(int m, int n);

	// 		 - ZEROSS FUNCTIONS
	void zeros();
	void zeros(int m, int n);


	//::M-03 - Arithmetic operations
	//	- ADD / SUBRACT / MULTIPLY
	void add(double v);
	void sub(double v);
	void mul(double v);
	void div(double v);

	//::M-04 - Utility functions
	//	- resize
	void resize(int m, int n);
	void add_column();
	void add_column(int j);
	void add_row();
	void add_row(int i);

	void delete_column(int j);
	void delete_row(int a);

	//::M-05 - minor
	//		 - MINOR FUNCTION
	CMatrix_ver2 minor_mk(const int r, const int c);
	void minor_mk2(const int r, const int c);
};


///////////////////////////////////
/* Global functions declarations */
///////////////////////////////////
CMatrix_ver2	operator+ (const CMatrix_ver2& A, 	const CMatrix_ver2& B);
CMatrix_ver2	operator+ (const CMatrix_ver2& A, 	const double B);
CMatrix_ver2	operator+ (const double B, 			const CMatrix_ver2& A);

CMatrix_ver2	operator- (const CMatrix_ver2& A, 	const CMatrix_ver2& B);
CMatrix_ver2	operator- (const CMatrix_ver2& A, 	const double B);
CMatrix_ver2	operator- (const double B, 			const CMatrix_ver2& A);
CMatrix_ver2	operator- (const CMatrix_ver2& A);

CMatrix_ver2	operator* (const CMatrix_ver2& A, 	const CMatrix_ver2& B);
CMatrix_ver2	operator* (const double B,			const CMatrix_ver2& A);

CMatrix_ver2	CMatrix_ZEROS(int i, int j);
CMatrix_ver2 	CMatrix_ONES(int i, int j);
CMatrix_ver2 	CMatrix_DIAG(int i);

CMatrix_ver2	CMatrix_INV(CMatrix_ver2& A);
double 			CMatrix_DET(CMatrix_ver2& A);


#endif /* CMATRIX_HPP_ */
