/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 28, 2015
 *      			Author: Minkwan Kim
 *
 * FaceGeometryBasic.hpp
 * 			-  
 *  
 */
#ifndef FACEGEOMETRYBASIC_HPP_
#define FACEGEOMETRYBASIC_HPP_

#include <iostream>
#include <vector>

using namespace std;


namespace OP2A{
namespace GRID{

class Node;

enum FaceType
{
	f_mixed				= 0,
	f_line				= 2,
	f_triangle			= 3,
	f_quadrilateral 	= 4
};



/*
 * Class for face basic geometry
 * @author Minkwan Kim
 * @version 1.0 27/5/2015
 */

class FaceGeoBasic
{
public:
	FaceType	type;

	int NN;
	vector<Node *> node_list;


	double S;
	vector<double> x;
	vector< vector<double > > n;

	unsigned int BC;
	double	dist_wall;
	double	n_dot_wall;

	FaceGeoBasic();
	explicit FaceGeoBasic(const int ND, const FaceType f_type);

	~FaceGeoBasic();

	void allocate(const int ND, const FaceType f_type);
	bool allocation();

private:
	bool is_allocated;
};

}
}



#endif /* FACEGEOMETRYBASIC_HPP_ */
