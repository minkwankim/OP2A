/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 27, 2015
 *      			Author: Minkwan Kim
 *
 * CellGeometryCart.hpp
 * 			-  
 *  
 */
#ifndef CELLGEOMETRYCART_HPP_
#define CELLGEOMETRYCART_HPP_




#include <vector>
using namespace std;


namespace OP2A{
namespace GRID{

class CELL;

enum CellTypeCart
{
	flow	= 0,
	surface	= 1,
	solid	= 2
};

class CellGeoCart
{
public:
	CellTypeCart	typecart;
	int 			cut_in_cell_type;

	int grid_level;
	bool need_to_refine;
	bool has_children;
	vector<CELL*>	children_list;

	CellGeoCart();
	~CellGeoCart();
};



}
}

#endif /* CELLGEOMETRYCART_HPP_ */
