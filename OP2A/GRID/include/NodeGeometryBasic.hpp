/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 28, 2015
 *      			Author: Minkwan Kim
 *
 * NodeGeometryBasic.hpp
 * 			-  
 *  
 */
#ifndef NODEGEOMETRYBASIC_HPP_
#define NODEGEOMETRYBASIC_HPP_


#include <vector>
using namespace std;


namespace OP2A{
namespace GRID{

class NodeGeoBasic
{
public:
	int	ID;
	vector<double> x;

	NodeGeoBasic():x(2, 0.0), ID(-1)
	{

	}

	~NodeGeoBasic(){	};


protected:
	explicit NodeGeoBasic(const int ND): x(ND, 0.0), ID(-1)
	{

	}

	void reassignBasic(const int ND)
	{
		x.resize(ND);
		for (int i = 0; i <= ND-1; i++)	x[i] = 0.0;
	}

};



}
}


#endif /* NODEGEOMETRYBASIC_HPP_ */
