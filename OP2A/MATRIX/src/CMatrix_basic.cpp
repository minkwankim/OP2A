/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Mar 23, 2015
 *      			Author: Minkwan Kim
 *
 * CMatrix_basic.cpp
 * 			-  
 *  
 */

#include "../include/CMatrix_basic.hpp"



void swap(double& a, double& b)
{
	double temp;

	temp	= a;
	a 		= b;
	b		= temp;
}


int	Mat_index_2D_to_1D(int i, int j, int M, int N)
{
	if (i <= 0)	Exception("[ERROR][matrix_MK]: i-index should be positive integer");
	if (j <= 0)	Exception("[ERROR][matrix_MK]: j-index should be positive integer");

	int k;
	k	= N*(i-1) + (j-1);

	return (k);
}
