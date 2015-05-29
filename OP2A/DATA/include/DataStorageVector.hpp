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

using namespace std;


namespace OP2A{
namespace Data{

class DataStorageVector
{
public:
	int numDataVector;
	vector<Data::DataStorage>		data;
	Common::Map1D<string, int>		dataMap;

	DataStorageVector();
	explicit DataStorageVector(const unsigned int size_data);
	explicit DataStorageVector(const unsigned int size_data, Common::Map1D<string, int> &data_map);


	void clean();
	void resize(unsigned int new_size);

	const Data::DataStorage& operator() (const unsigned int i);

	~DataStorageVector();

private:
	bool is_allocated;


};

}
}
#endif /* DATASTORAGEVECTOR_HPP_ */
