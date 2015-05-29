/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 4, 2015
 *      			Author: Minkwan Kim
 *
 * Cell.cpp
 * 			-  
 *  
 */



#include <vector>
using namespace std;

#include "../include/Cell.hpp"
#include "DATA/include/DataStorage.hpp"

namespace OP2A{
namespace GRID{


Cell::Cell(const int ND, const CellType c_type,	Data::DataStorageVector& data_template):geo(ND, c_type), data(data_template)
{

}


Cell::Cell(CellGeometry& Cellgeo_template, Data::DataStorageVector& data_template):geo(Cellgeo_template), data(data_template)
{

}



Cell::~Cell()
{

}

}
}
