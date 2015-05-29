/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 27, 2015
 *      			Author: Minkwan Kim
 *
 * CellGeometryCart.cpp
 * 			-  
 *  
 */



#include "../include/CellGeometryCart.hpp"

namespace OP2A{
namespace GRID{

CellGeoCart::CellGeoCart()
{
	typecart			= CellTypeCart::flow;
	cut_in_cell_type	= -1;

	grid_level	= 0;
	need_to_refine	= false;
	has_children	= false;
}

CellGeoCart::~CellGeoCart()
{

}



}
}
