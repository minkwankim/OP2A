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


Face::Face(const int ND, const FaceType f_type, bool extendedStencil,	Data::DataStorageVector<Data::DataStorage>& data_template)
:geo(ND, f_type, extendedStencil), data1D(data_template), allocated_geo(true), allocated_data1D(true), allocated_data2D(false)
{

}

Face::Face(FaceGeometry& geo_template,	Data::DataStorageVector<Data::DataStorage>& data_template)
:geo(geo_template), data1D(data_template), allocated_geo(true), allocated_data1D(true), allocated_data2D(false)
{

}


Face::Face(const int ND, const FaceType f_type, bool extendedStencil,	Data::DataStorageVector<Data::DataStorage>& data1D_template, Data::DataStorageVector<Data::DataStorage2D>& data2D_template)
:geo(ND, f_type, extendedStencil), data1D(data1D_template), data2D(data2D_template), allocated_geo(true), allocated_data1D(true), allocated_data2D(true)
{

}

Face::Face(FaceGeometry& geo_template, Data::DataStorageVector<Data::DataStorage>& data1D_template, Data::DataStorageVector<Data::DataStorage2D>& data2D_template)
:geo(geo_template), data1D(data1D_template), data2D(data2D_template), allocated_geo(true), allocated_data1D(true), allocated_data2D(true)
{

}




Face::~Face()
{

}


/*
 * internal functions
 */
void Face::allocate_geo_type(const int ND, const FaceType f_type, bool extendedStencil)
{
	geo.allocate(ND, f_type);
	if (extendedStencil == true) geo.apply_extended_sencil();

	allocated_geo = true;
}

void Face::allocate_geo_type(FaceGeometry& geo_template)
{
	geo	= geo_template;
	allocated_geo = true;
}

void Face::allocate_data1D_type(Data::DataStorageVector<Data::DataStorage>& data1D_template)
{
	data1D	= data1D_template;
	allocated_data1D = true;
}

void Face::allocate_data2D_type(Data::DataStorageVector<Data::DataStorage2D>& data2D_template)
{
	data2D	= data2D_template;
	allocated_data2D = true;
}


bool Face::is_allocated_geo()
{
	return allocated_geo;
}

bool Face::is_allocated_data1D()
{
	return allocated_data1D;
}

bool Face::is_allocated_data2D()
{
	return allocated_data2D;
}



}
}
