/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 11, 2015
 *      			Author: Minkwan Kim
 *
 * DataStorage2D.cpp
 * 			-  
 *  
 */



#include "../include/DataStorage2D.hpp"
#include "../include/Exception_DataStorageSize.hpp"


namespace OP2A{
namespace Data{



DataStorage2D::DataStorage2D():numData_I(1), numData_J(1),dataMap(1)
{

}

DataStorage2D::DataStorage2D(const unsigned data_size_I,  const unsigned data_size_J):
		numData_I(data_size_I), numData_J(data_size_J), data(data_size_I, data_size_J), dataMap(data_size_I*data_size_J)
{

}

DataStorage2D::DataStorage2D(string data_name, const unsigned data_size_I,  const unsigned data_size_J):
		name(data_name), numData_I(data_size_I), numData_J(data_size_J), data(data_size_I, data_size_J), dataMap(data_size_I*data_size_J)
{

}


DataStorage2D::DataStorage2D(string data_name, const unsigned data_size_I,  const unsigned data_size_J, Common::Map2D<string, string, int>& data_map):
		name(data_name), numData_I(data_size_I), numData_J(data_size_J), data(data_size_I, data_size_J), dataMap(data_map)
{

}


DataStorage2D::~DataStorage2D()
{

}



void DataStorage2D::resize(const unsigned int numdata_I, const unsigned int numdata_J)
{
	numData_I	= numdata_I;
	numData_J	= numdata_J;

	data.resize(numData_I, numData_J);
	dataMap.reserve(numData_I * numData_J);
}


void DataStorage2D::resize(const unsigned int numdata_I, const unsigned int numdata_J, Common::Map2D<string, string, int>& data_map)
{
	numData_I	= numdata_I;
	numData_J	= numdata_J;
	data.resize(numData_I, numData_J);

	dataMap = data_map;

	if (numData_I*numData_J != dataMap.size())
	{
		throw ExceptionDataStorageSize(FromHere(), "DataStroage2D size does not match with mapping data");
	}
}

void DataStorage2D::asgName(const string name_data)
{
	name	= name_data;
}


double& DataStorage2D::operator()(const unsigned i, const unsigned j)
{
	return data(i+1,j+1);
}

double& DataStorage2D::operator()(const string col, const string row)
{
	int index = dataMap.find(col, row);

	int i, j;
	i	= (index+1) 	% numData_I;
	j	= (index+1 - i) / numData_I + 1;

	return data(i,j);
}

void DataStorage2D::mapInsert(const string col, const string row, const unsigned i, const unsigned j)
{
	op_assert (i+1 <= numData_I);
	op_assert (j+1 <= numData_J);

	int	index = (i + numData_I*j);

	dataMap.insert(col, row, index);
}



}
}
