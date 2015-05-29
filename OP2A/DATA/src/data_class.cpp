/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 12, 2015
 *      			Author: Minkwan Kim
 *
 * data_class.cpp
 * 			-  
 *  
 */




#include <vector>

#include "../include/data_class.hpp"


using namespace std;


/*
 * Basic Data class
 */
DATA_BASIC_ver1::DATA_BASIC_ver1()
{
	for (int i = 0; i <= 3; i++)	NV[i]	= 0;
}

DATA_BASIC_ver1::~DATA_BASIC_ver1()
{

}

void DATA_BASIC_ver1::allocate_size()
{
	Q.resize(NV[0]);
	V.resize(NV[1]);
	W.resize(NV[2]);
	data.resize(NV[3]);
}


DATA_BASIC_ver2::DATA_BASIC_ver2()
{
	NV	= 0;
	ID	= 0;
}

DATA_BASIC_ver2::~DATA_BASIC_ver2()
{

}

void DATA_BASIC_ver2::allocate_size()
{
	data.resize(NV);
}


DATA_BASIC_2D::DATA_BASIC_2D()
{
	ID	= 0;
}

DATA_BASIC_2D::~DATA_BASIC_2D()
{

}

void DATA_BASIC_2D::allocate_size(int n, int m)
{
	data	= vector_2D(n, m, 0.0);
}





/*
 * Solution Data class
 */
SOL_CLASS_DATA::SOL_CLASS_DATA()
{

}

SOL_CLASS_DATA::~SOL_CLASS_DATA()
{
	int num	= data_ptr.size();

	//for (int n = 0; n <= num-1; n++)	data_ptr[n]	= NULL;

	for (int n = 0; n <= num-1; n++)	delete (data_ptr[n]);
	data_ptr.clear();
}

// Internal function
void SOL_CLASS_DATA::resize(unsigned int ncm, unsigned int nv)
{
	whereis.resize(ncm+1);

	for (int i = 0; i <= ncm-1; i++)
	{
		data_ptr.push_back(new DATA_BASIC());
		data_ptr[i]->NV	= nv;
		data_ptr[i]->allocate_size();
	}
}


SOL_CLASS_DATA_2D::SOL_CLASS_DATA_2D()
{

}

SOL_CLASS_DATA_2D::~SOL_CLASS_DATA_2D()
{
	int num	= data_ptr.size();

	for (int n = 0; n <= num-1; n++)	delete (data_ptr[n]);
	data_ptr.clear();
}

// Internal function
void SOL_CLASS_DATA_2D::resize(unsigned int ncm, unsigned int n, int m)
{
	whereis.resize(ncm+1);

	for (int i = 0; i <= ncm-1; i++)
	{
		data_ptr.push_back(new DATA_BASIC_2D());
		data_ptr[i]->allocate_size(n, m);
	}
}















/*
 * Solution class
 */
SOL_CLASS_BASIC::SOL_CLASS_BASIC()
{
	NDATA	= 0;
	NV		= 0;
};

SOL_CLASS_BASIC::~SOL_CLASS_BASIC()
{

}

void SOL_CLASS_BASIC::allocate_size(unsigned int ncm, unsigned int nv)
{
	NV		= nv;
	NDATA	= ncm;

	variable_names.resize(NV);
	resize(ncm, nv);
}


