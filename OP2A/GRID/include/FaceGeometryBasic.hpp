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

class FACE;
class NODE;

enum FaceType
{
	line			= 2,
	triangle		= 3,
	quadrilateral 	= 4
};

enum BCType
{
	interior		= 0,
	wall			= 1,
	inlet			= 2,
	outlet			= 3,
	freestream		= 4,
	symmetric		= 5,
	axis			= 6,
	anode			= 7,
	cathode			= 8,
	dielectricwall	= 9
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
	vector<NODE *> node_list;


	double S;
	vector<double> x;
	vector< vector<double > > n;

	unsigned int BC;
	double	dist_wall;
	double	n_dot_wall;

	~FaceGeoBasic();

protected:
	explicit FaceGeoBasic(const int ND, const FaceType f_type);
};

}
}



#endif /* FACEGEOMETRYBASIC_HPP_ */
