/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 4, 2015
 *      			Author: Minkwan Kim
 *
 * Cell.hpp
 * 			-  
 *  
 */
#ifndef CELL_HPP_
#define CELL_HPP_


#include <vector>
using namespace std;

#include "CellGeometry.hpp"

#include "DATA/include/DataStorage.hpp"
#include "DATA/include/DataStorageVector.hpp"

namespace OP2A{
namespace GRID{


/*
 * Class for Cell
 * @author Minkwan Kim
 * @version 1.0 27/5/2015
 */

#ifdef	USE_CART
#	define CellGeometry	CellGeometry2
#else
#	define CellGeometry	CellGeometry1
#endif


class Cell
{
public:
	CellGeometry				geo;
	Data::DataStorageVector		data;

	Cell()
	{

	};

	explicit	Cell(const int ND, const CellType c_type,	Data::DataStorageVector& data_template);
	explicit	Cell(CellGeometry& Cellgeo_template, Data::DataStorageVector& data_template);

	~Cell();
};






}
}


#endif /* CELL_HPP_ */
