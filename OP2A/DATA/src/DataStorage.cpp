/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 29, 2015
 *      			Author: Minkwan Kim
 *
 * DataStorage.cpp
 * 			-  
 *  
 */



#include "../include/DataStorage.hpp"
#include "../include/Exception_DataStorageSize.hpp"


namespace OP2A{
namespace Data{


DataStorage::DataStorage(string data_name, const unsigned data_size):name(data_name), numData(data_size), data(data_size, 0.0), dataMap(data_size)
{

}

DataStorage::DataStorage(string data_name, const unsigned data_size, Common::Map1D<string, int>	data_map):name(data_name), numData(data_size), data(data_size, 0.0), dataMap(data_map)
{
	if (numData != dataMap.size())
	{
		throw ExceptionDataStorageSize(FromHere(), "DataStroage size does not match with mapping data");
	}
}


DataStorage::~DataStorage()
{

}



}
}
