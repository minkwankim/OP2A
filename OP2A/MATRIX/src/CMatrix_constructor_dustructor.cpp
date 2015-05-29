/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Mar 23, 2015
 *      			Author: Minkwan Kim
 *
 * CMatrix_constructor_dustructor.cpp
 * 			-  
 *  
 */

#include "../include/CMatrix.hpp"



// A. Constructors
CMatrix_ver2::CMatrix_ver2()
{
	M 			= 0;
	N 			= 0;
	elements	= NULL;
}

CMatrix_ver2::CMatrix_ver2(int m, int n)
{
	M 			= m;
	N 			= n;
	elements	= new double [m*n];
}


CMatrix_ver2::CMatrix_ver2(vector<vector<double> > &other)
{
	M 			= other.size();
	N 			= other[0].size();
	elements	= new double [M*N];

	int k = 0;
	for (int i = 0; i <= M-1; i++)
	{
		for (int j = 0; j <= N-1; j++)
		{
			elements[k]	= other[i][j];
			k++;
		}
	}
}


// B. Deconstructors
CMatrix_ver2::~CMatrix_ver2()
{
	if (M*N > 0)
	{
		delete [] elements;
		N	= 0;
		M	= 0;
	}

	elements = NULL;
}
