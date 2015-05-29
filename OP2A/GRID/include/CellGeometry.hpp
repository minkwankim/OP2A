/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 27, 2015
 *      			Author: Minkwan Kim
 *
 * CellGeometry.hpp
 * 			-  
 *  
 */
#ifndef CELLGEOMETRY_HPP_
#define CELLGEOMETRY_HPP_


#include <vector>
using namespace std;

#include "CellGeometryBasic.hpp"
#include "CellGeometryCart.hpp"


namespace OP2A{
namespace GRID{

class	CellGeometry1 :public CellGeoBasic
{
public:
	explicit CellGeometry1(const int ND, const CellType c_type) : CellGeoBasic(ND, c_type)
	{

	}

	~CellGeometry1()
	{

	}
};


class	CellGeometry2 :public CellGeoBasic, public CellGeoCart
{
public:
	explicit CellGeometry2(const int ND, const CellType c_type) : CellGeoBasic(ND, c_type)
	{

	}

	~CellGeometry2()
	{

	}
};



}
}
#endif /* CELLGEOMETRY_HPP_ */
