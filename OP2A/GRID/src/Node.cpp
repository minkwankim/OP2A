/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 1, 2015
 *      			Author: Minkwan Kim
 *
 * Node.cpp
 * 			-  
 *  
 */



#include <vector>
using namespace std;

#include "../include/Node.hpp"

namespace OP2A{
namespace GRID{


Node::Node()
{

}

Node::Node(const int ND, const unsigned int nsc, Data::DataStorage	&data_node):geo(ND, nsc), data(data_node)
{

}

Node::Node(Node& node_sample):geo(node_sample.geo), data(node_sample.data)
{

}


void Node::assign(const unsigned int ND, const unsigned int nsc, Data::DataStorage &data_node)
{
	geo.assign(ND, nsc);
	data = data_node;
}


Node::~Node()
{

}





}
}
