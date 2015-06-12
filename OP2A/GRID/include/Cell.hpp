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
#include "DATA/include/DataStorage2D.hpp"
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
	CellGeometry									geo;
	Data::DataStorageVector<Data::DataStorage>		data1D;
	Data::DataStorageVector<Data::DataStorage2D>	data2D;


	Cell():allocated_geo(false), allocated_data1D(false), allocated_data2D(false)
	{

	};

	explicit	Cell(const int ND, const CellType c_type,	Data::DataStorageVector<Data::DataStorage>& data_template);
	explicit	Cell(CellGeometry& Cellgeo_template, Data::DataStorageVector<Data::DataStorage>& data_template);

	explicit	Cell(const int ND, const CellType c_type,	Data::DataStorageVector<Data::DataStorage>& data_template, Data::DataStorageVector<Data::DataStorage2D>& data2D_template);
	explicit	Cell(CellGeometry& Cellgeo_template, Data::DataStorageVector<Data::DataStorage>& data_template, Data::DataStorageVector<Data::DataStorage2D>& data2D_template);

	~Cell();


public:
	void allocate_geo_type(const int ND, const CellType c_type);
	void allocate_geo_type(CellGeometry& geo_template);
	void allocate_data1D_type(Data::DataStorageVector<Data::DataStorage>& data1D_template);
	void allocate_data2D_type(Data::DataStorageVector<Data::DataStorage2D>& data2D_template);

	bool is_allocated_geo();
	bool is_allocated_data1D();
	bool is_allocated_data2D();


private:
	bool allocated_geo;
	bool allocated_data1D;
	bool allocated_data2D;
};






}
}


#endif /* CELL_HPP_ */
