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

class Face;
class Node;

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
	vector<Node *> node_list;

	int	NF;
	vector<Face *> face_list;

	double S;
	vector<double> x;

	double	characteristic_length;			// characteristic_length
	double 	dist_wall;						// distance to wall


	CellGeoBasic();
	explicit CellGeoBasic(const int ND, const CellType c_type);

	virtual ~CellGeoBasic() {	}

public:
	void allocate(const int ND, const CellType c_type);
	bool is_allocated();


private:
	bool allocated;

};





}
}





#endif /* CELLGEOMETRYBASIC_HPP_ */
