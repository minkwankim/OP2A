/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 28, 2015
 *      			Author: Minkwan Kim
 *
 * NodeGeometry.hpp
 * 			-  
 *  
 */
#ifndef NODEGEOMETRY_HPP_
#define NODEGEOMETRY_HPP_

#include <vector>
using namespace std;

#include "NodeGeometryBasic.hpp"
#include "NodeGeometryAdd.hpp"

namespace OP2A{
namespace GRID{

class NodeGeometry : public NodeGeoBasic , public NodeGeoAdd
{
public:
	NodeGeometry()
	{

	}

	explicit NodeGeometry(const int ND) : NodeGeoBasic(ND)
	{

	}

	explicit NodeGeometry(const int ND, unsigned int const nsc) : NodeGeoBasic(ND), NodeGeoAdd(nsc)
	{

	}

	void assign(const int ND, unsigned int const nsc)
	{
		reassignBasic(ND);
		reassignAdd(nsc);
	}



	~NodeGeometry()	{	}
};



}
}
#endif /* NODEGEOMETRY_HPP_ */
