/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 2, 2015
 *      			Author: Minkwan Kim
 *
 * DataStorageVector.hpp
 * 			-  
 *  
 */
#ifndef DATASTORAGEVECTOR_HPP_
#define DATASTORAGEVECTOR_HPP_


#include <vector>

#include "DATA/include/DataStorage.hpp"
#include "DATA/include/Exception_DataStorageVector.hpp"

using namespace std;


namespace OP2A{
namespace Data{


template <typename TYPE>
class DataStorageVector
{
public:
	int numDataVector;
	vector<TYPE>					data;
	Common::Map1D<string, int>		dataMap;

	DataStorageVector():numDataVector(0), dataMap(1), is_allocated(false), is_mapped(false)
	{

	}

	explicit DataStorageVector(const unsigned int size_data):numDataVector(size_data), dataMap(size_data), is_mapped(false)
	{
		data.resize(numDataVector);
		is_allocated	= true;
	}


	explicit DataStorageVector(const unsigned int size_data, Common::Map1D<string, int> &data_map):numDataVector(size_data), dataMap(data_map)
	{
		data.resize(numDataVector);
		is_allocated	= true;
		is_mapped		= true;

		if (numDataVector != dataMap.size())	throw ExceptionDataStorageSize(FromHere(), "DataStroageVector size does not match with mapping data");
	}



	void clean()
	{
		if (is_allocated == true)
		{
			data.clear();

			is_allocated = false;
			is_mapped	 = false;

			numDataVector	= 0;
			dataMap.clear();
		}
		else
		{
			numDataVector	= 0;

			is_mapped	= false;
			dataMap.clear();
		}
	}



	void resize(unsigned int new_size)
	{
		data.resize(new_size);
		numDataVector = new_size;
		dataMap.reserve(numDataVector);

		is_mapped	= false;
	}


	void mapping()
	{
		dataMap.clear();
		dataMap.reserve(numDataVector);

		for (int i = 0; i <= numDataVector-1; i++)	dataMap.insert(data[i].name, i);

		is_mapped	= true;
	}


	TYPE& operator() (const unsigned int i)
	{
		if (i >= numDataVector)
		{
			throw ExceptionDataStorageVector(FromHere(), "Exceed the data size of DataStorageVector");
		}

		return data[i];
	}




	TYPE& operator()(const string var_main)
	{
		int index_main	= dataMap.find(var_main);


		return data[index_main];
	}


	double& operator()(const string var_main, const string var1)
	{
		int index_main	= dataMap.find(var_main);


		return data[index_main](var1);
	}

	double& operator()(const string var_main, const string var1, const string var2)
	{
		int index_main	= dataMap.find(var_main);


		return data[index_main](var1, var2);
	}











	~DataStorageVector()
	{

	}

private:
	bool is_allocated;
	bool is_mapped;


};








}
}
#endif /* DATASTORAGEVECTOR_HPP_ */
