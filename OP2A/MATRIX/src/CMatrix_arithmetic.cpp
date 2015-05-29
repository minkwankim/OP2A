/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Mar 24, 2015
 *      			Author: Minkwan Kim
 *
 * CMatrix_arithmetic.cpp
 * 			-  
 *  
 */



#include "../include/CMatrix.hpp"


void CMatrix_ver2::add(double v)
{
	if (N*M == 0)	Exception("[ERROR][CMatrix]: Size is not defined yet.");

	for (int i = 0; i <= M*N-1; i++)	elements[i] += v;
}


void CMatrix_ver2::sub(double v)
{
	if (N*M == 0)	Exception("[ERROR][CMatrix]: Size is not defined yet.");

	for (int i = 0; i <= M*N-1; i++)	elements[i] -= v;
}


void CMatrix_ver2::mul(double v)
{
	if (N*M == 0)	Exception("[ERROR][CMatrix]: Size is not defined yet.");

	for (int i = 0; i <= M*N-1; i++)	elements[i] *= v;
}

void CMatrix_ver2::div(double v)
{
	if (N*M == 0)	Exception("[ERROR][CMatrix]: Size is not defined yet.");

	for (int i = 0; i <= M*N-1; i++)	elements[i] /= v;
}
