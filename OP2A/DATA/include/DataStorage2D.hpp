/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 11, 2015
 *      			Author: Minkwan Kim
 *
 * DataStorage2D.hpp
 * 			-  
 *  
 */
#ifndef DATASTORAGE2D_HPP_
#define DATASTORAGE2D_HPP_



#include <vector>

#include "Common/include/Map1D.hpp"
#include "Common/include/Map2D.hpp"
#include "Common/include/Map3D.hpp"
#include "Common/include/Vector2D.hpp"

using namespace std;


namespace OP2A{
namespace Data{

class DataStorage2D
{
public:
	string						name;
	unsigned int				numData_I;
	unsigned int				numData_J;
	Common::Vector2D<double>	data;


	DataStorage2D();
	explicit DataStorage2D(const unsigned data_size_I,  const unsigned data_size_J);
	explicit DataStorage2D(string data_name, const unsigned data_size_I,  const unsigned data_size_J);
	explicit DataStorage2D(string data_name, const unsigned data_size_I,  const unsigned data_size_J, Common::Map2D<string, string, int>& data_map);


	void resize(const unsigned int numdata_I, const unsigned int numdata_J);
	void resize(const unsigned int numdata_I, const unsigned int numdata_J, Common::Map2D<string, string, int>& data_map);
	void asgName(const string name_data);

	double& operator()(const unsigned i, const unsigned j);
	double& operator()(const string col, const string row);
	void mapInsert(const string col, const string row, const unsigned i, const unsigned j);

	~DataStorage2D();


private:
	Common::Map2D<string, string, int>	dataMap;
};

}
}

#endif /* DATASTORAGE2D_HPP_ */
