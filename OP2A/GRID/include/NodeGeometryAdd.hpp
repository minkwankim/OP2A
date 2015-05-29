/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 28, 2015
 *      			Author: Minkwan Kim
 *
 * NodeGeometryAdd.hpp
 * 			-  
 *  
 */
#ifndef NODEGEOMETRYADD_HPP_
#define NODEGEOMETRYADD_HPP_



#include <vector>
using namespace std;


namespace OP2A{
namespace GRID{

class CELL;


class NodeGeoAdd
{
public:
	unsigned int	NSC;
	vector<CELL *>	CellList;
	vector<double> 	WeightShairedCell;

	NodeGeoAdd()
	{
		NSC = 0;
	}

	NodeGeoAdd(unsigned int const nsc) : NSC(nsc), CellList(vector<CELL*> [nsc]), WeightShairedCell(vector<double> [nsc, 0.0])
	{

	}

	virtual ~NodeGeoAdd();
};



}
}



#endif /* NODEGEOMETRYADD_HPP_ */
