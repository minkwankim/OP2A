/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 1, 2015
 *      			Author: Minkwan Kim
 *
 * Face.cpp
 * 			-  
 *  
 */




#include <vector>
using namespace std;

#include "../include/Face.hpp"
#include "DATA/include/DataStorage.hpp"

namespace OP2A{
namespace GRID{


Face::Face(const int ND, const FaceType f_type, bool extendedStencil,	Data::DataStorageVector& data_template):geo(ND, f_type, extendedStencil), data(data_template)
{

}

Face::Face(FaceGeometry& geo_template,	Data::DataStorageVector& data_template):geo(geo_template), data(data_template)
{

}



Face::~Face()
{

}

}
}
