/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 1, 2015
 *      			Author: Minkwan Kim
 *
 * Face.hpp
 * 			-  
 *  
 */
#ifndef FACE_HPP_
#define FACE_HPP_



#include <vector>
using namespace std;

#include "FaceGeometry.hpp"

#include "DATA/include/DataStorage.hpp"
#include "DATA/include/DataStorageVector.hpp"

namespace OP2A{
namespace GRID{


/*
 * Class for Face
 * @author Minkwan Kim
 * @version 1.0 27/5/2015
 */

class Face
{
public:
	FaceGeometry				geo;
	Data::DataStorageVector		data;

	Face() {	};

	explicit	Face(const int ND, const FaceType f_type, bool extendedStencil,	Data::DataStorageVector& data_template);
	explicit	Face(FaceGeometry& geo_template,	Data::DataStorageVector& data_template);

	~Face();
};








}
}
#endif /* FACE_HPP_ */
