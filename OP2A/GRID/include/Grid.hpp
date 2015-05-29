/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 29, 2015
 *      			Author: Minkwan Kim
 *
 * Grid.hpp
 * 			-  
 *  
 */
#ifndef GRID_HPP_
#define GRID_HPP_


#include "Node.hpp"
#include "Face.hpp"
#include "Cell.hpp"

#include "Common/include/Map1D.hpp"
#include "Common/include/Map2D.hpp"
#include "Common/include/Map3D.hpp"


namespace OP2A{
namespace GRID{


class GridInfO {
public:
    const unsigned int NCM;
    const unsigned int NFM;
    const unsigned int NCM;
    const unsigned int NGM;
    

public:
  <#member functions#>
};

class Grid
{
public:
	unsigned int NNM;
	unsigned int NFM;
	unsigned int NCM;
	unsigned int NGM;

	vector<Node>	nodes;
	vector<Face>	faces;
	vector<Cell>	cells;
	vector<Cell>	cells_ghost;

public:	// Constructor and desctructor
	Grid();
	~Grid();


public:	// Functions


};







}
}



#endif /* GRID_HPP_ */
