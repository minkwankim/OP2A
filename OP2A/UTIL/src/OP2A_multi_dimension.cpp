/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 17, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_multi_dimension.cpp
 * 			-  
 *  
 */




#include "../include/OP2A_multi_dimension.hpp"

vector< vector <double*> > 				NEW_vector_2D_ptr(int I, int J)
{
	vector< vector <double*> > vector_2D_ptr(I, vector<double*>(J));
	return(vector_2D_ptr);
}

vector< vector< vector <double*> > > 	NEW_vector_3D_ptr(int I, int J, int K)
{
	vector< vector< vector <double*> > > vector_3D_ptr(I, vector<vector<double*> >(J, vector<double*>(K)));
	return(vector_3D_ptr);
}




void DELETE_vector_2D_ptr(vector< vector <double*> > vector2D_ptr, int I, int J)
{
	for (int i = 0; i <= I-1; i++)
		for (int j = 0; j <= J-1; j++)
			delete vector2D_ptr[i][j];
}


void DELETE_vector_3D_ptr(vector< vector< vector <double*> > > vector3D_ptr, int I, int J, int K)
{
	for (int i = 0; i <= I-1; i++)
		for (int j = 0; j <= J-1; j++)
			for (int k = 0; k <= K-1; k++)
				delete vector3D_ptr[i][j][k];
}
