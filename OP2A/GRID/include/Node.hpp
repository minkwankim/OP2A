/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 29, 2015
 *      			Author: Minkwan Kim
 *
 * Node.hpp
 * 			-  
 *  
 */
#ifndef NODE_HPP_
#define NODE_HPP_



#include <vector>
using namespace std;

#include "NodeGeometry.hpp"
#include "DATA/include/DataStorage.hpp"

namespace OP2A{
namespace GRID{



/*
 * Class for Node
 * @author Minkwan Kim
 * @version 1.0 27/5/2015
 */

class Node
{
public:
	NodeGeometry			geo;
	Data::DataStorage		data;


	Node();
	explicit	Node(const int ND, const unsigned int nsc, Data::DataStorage	&data_node);
	explicit	Node(Node& node_sample);


	~Node();

	void assign(const unsigned int ND, const unsigned int nsc, Data::DataStorage &data_node);
};



}
}

#endif /* NODE_HPP_ */
