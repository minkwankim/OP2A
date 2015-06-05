/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 27, 2015
 *      			Author: Minkwan Kim
 *
 * CellGeometryBasic.cpp
 * 			-  
 *  
 */



#include "../include/CellGeometryBasic.hpp"
#include "../include/Exception_AllocationCell.hpp"


namespace OP2A{
namespace GRID{

CellGeoBasic::CellGeoBasic():ID(-1), type(ghost), NN(0), NF(0), S(0), characteristic_length(0), dist_wall(0), allocated(false), BC(0)
{

}


CellGeoBasic::CellGeoBasic(const int ND, const CellType c_type): ID(-1), type(c_type), x(ND, 0.0), allocated(true), BC(0)
{

	switch(type)
	{
	case CellType::ghost:
		NF	= 1;
		NN	= 4;
		break;

	case CellType::triangle:
		NF	= 3;
		NN	= 3;
		break;

	case CellType::tetrahedron:
		NF	= 4;
		NN	= 4;
		break;

	case CellType::quadrilateral:
		NF	= 4;
		NN	= 4;
		break;

	case CellType::hexahedron:
		NF	= 6;
		NN	= 8;
		break;

	case CellType::pryramid:
		NF	= 5;
		NN	= 5;
		break;

	case CellType::wedge:
		NF	= 5;
		NN	= 6;
		break;
	}

	node_list.resize(NN);
	face_list.resize(NF);

	S						= 0.0;
	characteristic_length	= 0.0;
	dist_wall				= 0.0;
}



void CellGeoBasic::allocate(const int ND, const CellType c_type)
{
	type = c_type;
	x.resize(ND, 0.0);
	allocated = true;

	switch(type)
	{
	case CellType::ghost:
		NF	= 1;
		NN	= 4;
		break;

	case CellType::triangle:
		NF	= 3;
		NN	= 3;
		break;

	case CellType::tetrahedron:
		NF	= 4;
		NN	= 4;
		break;

	case CellType::quadrilateral:
		NF	= 4;
		NN	= 4;
		break;

	case CellType::hexahedron:
		NF	= 6;
		NN	= 8;
		break;

	case CellType::pryramid:
		NF	= 5;
		NN	= 5;
		break;

	case CellType::wedge:
		NF	= 5;
		NN	= 6;
		break;
	}

	node_list.resize(NN);
	face_list.resize(NF);

	S						= 0.0;
	characteristic_length	= 0.0;
	dist_wall				= 0.0;
}

bool CellGeoBasic::is_allocated()
{
	return (allocated);
}



}
}
