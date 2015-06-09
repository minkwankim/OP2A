/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 9, 2015
 *      			Author: Minkwan Kim
 *
 * MultiDimention.hpp
 * 			-  
 *  
 */
#ifndef MULTIDIMENTION_HPP_
#define MULTIDIMENTION_HPP_


#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <cmath>
#include <memory>
#include "limits"
#include <vector>


#ifdef MPIF
	#include <mpi.h>
#endif

using namespace std;

namespace OP2A{
namespace Common{


/*
 * ================================
 * 	2D Vector
 * ================================
 */
template <class T>
vector<vector<T> > vector_2D(int n, int m, T value)
{
	vector<vector<T> > myVector(n, vector<T>(m, value));
	return myVector;
}

template <class T>
vector<vector<T> > vector_2D(int n, int m)
{
	vector<vector<T> > myVector(n, vector<T>(m));
	return myVector;
}


/*
 * ================================
 * 	3D Vector
 * ================================
 */
template <class T>
vector <vector<vector<T> > > vector_3D(int n, int m, int k, T value)
{
	vector <vector<vector<T> > > myVector(n, vector<vector<T> >(m, vector<T>(k, value)));

	return myVector;
}


template <class T>
vector <vector<vector<T> > > vector_3D(int n, int m, int k)
{
	vector <vector<vector<T> > > myVector(n, vector<vector<T> >(m, vector<T>(k)));

	return myVector;
}


/*
 * ================================
 * 	Resize Pointer and Array
 * ================================
 */
template <class T>
static void resize_class_pointer(T*& class_pointer, int N_previous, int N_new)
{
	if (N_previous < N_new)
	{
		T* ret	= new T[N_new];
		memcpy(ret, class_pointer, sizeof(class_pointer[0])*N_previous);
		delete[] class_pointer;
		class_pointer = ret;
	}
}

// Resize Array
template <class T>
void resize_array(T *cell_tree_c, int num_Tree_c, T **cell_tree_u, int num_Tree_u)
{
	int 		i_index;
	T	*cell_tree_temp;	// Temporary Tree cell array

	// Assign Pointer of new array
	*cell_tree_u	= new T [num_Tree_u+1];
	cell_tree_temp	= *cell_tree_u;						// Indicate pointer of updated cell array

	// Assign Previous cell data
	for (i_index = 1; i_index <= num_Tree_c; i_index++)
	{
		cell_tree_temp[i_index]	= cell_tree_c[i_index];
	}
}




/*
 * ================================
 * 	Multi-dimenstion Pointer and Array
 * ================================
 */

template <class T>
void DELETE_VEC(T* vec1d)
{
	vec1d	= NULL;
	delete[] vec1d;
}


template <class T>
T** NEW_MAT(int n, int m)
{
	T** buf2d;

	buf2d = new T* [n];

	for(int i = 0; i <= n-1; i++)
	{
		buf2d[i] = new T [m];
	}

	return buf2d;
}

template <class T>
void DELETE_MAT(T** buf2d, int n, int m)
{
	for(int i = 0; i <= n-1; i++)
	{
		delete[] buf2d[i];
	}

	delete[] buf2d;
	buf2d = NULL;
}


template <class T>
T*** NEW_MAT(int l, int m, int n)
{
	T*** buf3d;

	buf3d = new T** [l];

	for(int i = 0; i < l; i++)
	{
		buf3d[i] = new T* [m];
		for (int j = 0; j < m; j++)
		{
			buf3d[i][j] = new T [n];
		}
	}

	return buf3d;
}

template <class T>
void DELETE_MAT(T*** buf3d, int l, int m, int n)
{
	for(int i = 0; i <= l-1; i++)
	{
		for (int j = 0; j <= m-1; j++)
		{
			delete[] buf3d[i][j];
		}

		delete[] buf3d[i];
	}

	delete[] buf3d;
	buf3d = NULL;
}


template <class T>
void remove_element_in_vector(T value, vector<T> &data)
{
	int N	= data[0];
	for (int i = 1; i <= N; i++)
	{
		if (data[i] == value)
		{
			for (int j = i; j <= N-1; j++)	data[j]	= data[j+1];
			data[0]	= data[0] - 1;
		}
	}
}






/*
 * 2D/3D vector of pointers
 */
vector< vector <double*> > 				NEW_vector_2D_ptr(int I, int J);
vector< vector< vector <double*> > > 	NEW_vector_3D_ptr(int I, int J, int K);
void DELETE_vector_2D_ptr(vector< vector <double*> > vector2D_ptr, int I, int J);
void DELETE_vector_3D_ptr(vector< vector< vector <double*> > > vector3D_ptr, int I, int J, int K);



}
}


#endif /* MULTIDIMENTION_HPP_ */
