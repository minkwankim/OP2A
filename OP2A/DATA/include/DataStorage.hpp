/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 29, 2015
 *      			Author: Minkwan Kim
 *
 * DataStorage.hpp
 * 			-  
 *  
 */
#ifndef DATASTORAGE_HPP_
#define DATASTORAGE_HPP_


#include <vector>

#include "Common/include/Map1D.hpp"
#include "Common/include/Map2D.hpp"
#include "Common/include/Map3D.hpp"
#include "Common/include/Vector2D.hpp"

using namespace std;


namespace OP2A{
namespace Data{

class DataStorage
{
public:
	string				name;
	unsigned int		numData;
	vector<double>		data;

	Common::Map1D<string, int>	dataMap;

	DataStorage();
	explicit DataStorage(string data_name, const unsigned data_size);
	explicit DataStorage(string data_name, const unsigned data_size, Common::Map1D<string, int>& data_map);
	//explicit DataStorage(DataStorage &data_template);


	void resize(const unsigned int numdata);
	void resize(const unsigned int numdata, Common::Map1D<string, int>&	data_map);

	void asgName(const string name_data);

	~DataStorage();
};



}
}


#endif /* DATASTORAGE_HPP_ */
