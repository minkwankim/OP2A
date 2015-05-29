/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 27, 2015
 *      			Author: Minkwan Kim
 *
 * CellGeometryBasic.hpp
 * 			-  
 *  
 */
#ifndef CELLGEOMETRYBASIC_HPP_
#define CELLGEOMETRYBASIC_HPP_


#include <iostream>
#include <vector>

using namespace std;


namespace OP2A{
namespace GRID{

class FACE;
class NODE;

enum CellType
{
	ghost			= 0,
	triangle		= 1,
	tetrahedron		= 2,
	quadrilateral 	= 3,
	hexahedron		= 4,
	pryramid		= 5,
	wedge			= 6
};



/*
 * Class for cell basic geometry
 * @author Minkwan Kim
 * @version 1.0 27/5/2015
 */

class CellGeoBasic
{
public:
	CellType	type;

	int NN;
	vector<NODE *> node_list;

	int	NF;
	vector<FACE *> face_list;

	double S;
	vector<double> x;

	double	characteristic_length;			// characteristic_length
	double 	dist_wall;						// distance to wall


	virtual ~CellGeoBasic() {	}

protected:
	explicit CellGeoBasic(const int ND, const CellType c_type);

};





}
}





#endif /* CELLGEOMETRYBASIC_HPP_ */
