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



#include "../include/DataStorageVector.hpp"
#include "../include/Exception_DataStorageSize.hpp"
#include "../include/Exception_DataStorageVector.hpp"


namespace OP2A{
namespace Data{



DataStorageVector::DataStorageVector():numDataVector(0), dataMap(1), is_allocated(false), is_mapped(false)
{

}


DataStorageVector::DataStorageVector(const unsigned int size_data):numDataVector(size_data), dataMap(size_data), is_mapped(false)
{
	//data			= new Data::DataStorage	[numDataVector];
	data.resize(numDataVector);
	is_allocated	= true;
}


DataStorageVector::DataStorageVector(const unsigned int size_data, Common::Map1D<string, int> &data_map):numDataVector(size_data), dataMap(data_map)
{
	//data			= new Data::DataStorage	[numDataVector];
	data.resize(numDataVector);
	is_allocated	= true;
	is_mapped		= true;

	if (numDataVector != dataMap.size())	throw ExceptionDataStorageSize(FromHere(), "DataStroageVector size does not match with mapping data");
}


void DataStorageVector::clean()
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


void DataStorageVector::resize(unsigned int new_size)
{
	data.resize(new_size);
	numDataVector = new_size;
	dataMap.reserve(numDataVector);

	is_mapped	= false;
}




const Data::DataStorage &DataStorageVector::operator() (const unsigned int i)
{
	if (i >= numDataVector)
	{
		throw ExceptionDataStorageVector(FromHere(), "Exceed the data size of DataStorageVector");
	}

	return data[i];
}




DataStorageVector::~DataStorageVector()
{
	/*
	if (is_allocated == true)
	{
		delete [] data;
		data = NULL;

		is_allocated = false;
	}
	*/
}







/*
 * Internal Functions
 */
void DataStorageVector::mapping()
{
	dataMap.clear();
	dataMap.reserve(numDataVector);

	for (int i = 0; i <= numDataVector-1; i++)	dataMap.insert(data[i].name, i);

	is_mapped	= true;
}










}
}
