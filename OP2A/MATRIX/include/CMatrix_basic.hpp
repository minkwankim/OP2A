/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Mar 23, 2015
 *      			Author: Minkwan Kim
 *
 * CMatrix_basic.hpp
 * 			-  
 *  
 */

#ifndef CMATRIX_BASIC_HPP_
#define CMATRIX_BASIC_HPP_

#include "mkl.h"
#include <cmath>
#include <sstream>
#include <iostream>
#include <vector>
#include "../../UTIL/include/OP2A_utilities.hpp"


using namespace std;

/* EXCEPTION AND ERROR HANDLING */
class Exception{
	public:
	const char* msg;

	Exception(const char* arg) : msg(arg)
	{
		cerr << msg << endl;
		exit(0);
	};
};

void swap(double& a, double& b);
int	Mat_index_2D_to_1D(int i, int j, int M, int N);




#endif /* CMATRIX_BASIC_HPP_ */
