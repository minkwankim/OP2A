/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 3, 2015
 *      			Author: Minkwan Kim
 *
 * VectorMK.hpp
 * 			-  
 *  
 */
#ifndef VECTOR1D_HPP_
#define VECTOR1D_HPP_

#include <vector>
#include "Common/include/Exception_MemoryAllocator.hpp"
#include "GRID/include/Face.hpp"


namespace OP2A{
namespace Common{


template <typename TYPE>
class Vector1D
{
public:
	Vector1D()
	{
		num_size		= 0;
		is_allocated 	= false;
		data			= NULL;
	}

	explicit Vector1D(unsigned int i_size): num_size(i_size), is_allocated(true)
	{
		data	= new TYPE [i_size];
	}

	explicit Vector1D(unsigned int i_size, TYPE& data_sample): num_size(i_size), is_allocated(true)
	{
		data	= new TYPE [i_size];
		for (int i = 0; i <= num_size-1; i++)	data[i]	= data_sample;
	}

	~Vector1D()
	{
		if (is_allocated == true)	delete [] data;
		data = NULL;
	}



	void resize(unsigned int new_size)
	{
		TYPE* data_temp	=	new TYPE [new_size];

		for (int index = 0; index < num_size; index++)
		{
			data_temp[index]	= data[index];
		}

		num_size = new_size;

		delete [] data;
		data = data_temp;
	}

	unsigned int size()
	{
		return (num_size);
	}



	void push()
	{
		resize(num_size+1);
	}


	void push(TYPE& data_sample)
	{
		resize(num_size+1);
		data[num_size-1]	= data_sample;
	}

	void del(unsigned int i)
	{
		TYPE* data_temp	=	new TYPE [num_size-1];


		for (int index = 0; 	index < i; 			index++)	data_temp[index]	= data[index];
		for (int index = i+1; 	index < num_size;	index++)	data_temp[index-1]	= data[index];

		num_size--;

		delete [] data;
		data = data_temp;
	}

	void del()
	{
		TYPE* data_temp	=	new TYPE [num_size-1];


		for (int index = 0; 	index < num_size-1; 			index++)	data_temp[index]	= data[index];

		num_size--;

		delete [] data;
		data = data_temp;
	}




	const TYPE& operator[] (const unsigned int i)
	{
		if (i >= num_size)
		{
			throw ExceptionMemoryAllocator(FromHere(), "Exceed the data size of DataStorageVector");
		}

		return data[i];
	}






private:
	int num_size;
	TYPE* data;
	bool is_allocated;
};


}
}
#endif /* VECTORMK_HPP_ */
