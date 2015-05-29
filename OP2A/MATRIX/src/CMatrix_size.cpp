/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Mar 23, 2015
 *      			Author: Minkwan Kim
 *
 * CMatrix_size.cpp
 * 			-  
 *  
 */


#include "../include/CMatrix.hpp"


void CMatrix_ver2::size(int m, int n)
{
	if ( m <= 0 || n <= 0)	throw Exception("[ERROR][CMatrix]: Size of matrix should be unsigned integer!");

	M			= m;
	N 			= n;
	elements	= new double [M*N];
}

