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


Cell::Cell(const int ND, const CellType c_type,	Data::DataStorageVector<Data::DataStorage>& data_template)
:geo(ND, c_type), data1D(data_template), allocated_geo(true), allocated_data1D(true), allocated_data2D(false)
{

}


Cell::Cell(CellGeometry& Cellgeo_template, Data::DataStorageVector<Data::DataStorage>& data_template)
:geo(Cellgeo_template), data1D(data_template), allocated_geo(true), allocated_data1D(true), allocated_data2D(false)
{

}

Cell::Cell(const int ND, const CellType c_type,	Data::DataStorageVector<Data::DataStorage>& data1D_template,	Data::DataStorageVector<Data::DataStorage2D>& data2D_template)
:geo(ND, c_type), data1D(data1D_template), data2D(data2D_template), allocated_geo(true), allocated_data1D(true), allocated_data2D(true)
{

}

Cell::Cell(CellGeometry& Cellgeo_template, Data::DataStorageVector<Data::DataStorage>& data_template, 	Data::DataStorageVector<Data::DataStorage2D>& data2D_template)
:geo(Cellgeo_template), data1D(data_template), data2D(data2D_template), allocated_geo(true), allocated_data1D(true), allocated_data2D(true)
{

}



Cell::~Cell()
{

}






void Cell::allocate_geo_type(const int ND, const CellType c_type)
{
	geo.allocate(ND, c_type);
	allocated_geo = true;
}

void Cell::allocate_geo_type(CellGeometry& geo_template)
{
	geo	= geo_template;
	allocated_geo = true;
}

void Cell::allocate_data1D_type(Data::DataStorageVector<Data::DataStorage>& data1D_template)
{
	data1D	= data1D_template;
	allocated_data1D = true;
}

void Cell::allocate_data2D_type(Data::DataStorageVector<Data::DataStorage2D>& data2D_template)
{
	data2D	= data2D_template;
	allocated_data2D = true;
}


bool Cell::is_allocated_geo()
{
	return allocated_geo;
}

bool Cell::is_allocated_data1D()
{
	return allocated_data1D;
}

bool Cell::is_allocated_data2D()
{
	return allocated_data2D;
}


}
}
