/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 29, 2015
 *      			Author: Minkwan Kim
 *
 * FaceGeometry.hpp
 * 			-  
 *  
 */
#ifndef FACEGEOMETRY_HPP_
#define FACEGEOMETRY_HPP_

#include "FaceGeometryBasic.hpp"
#include "FaceStencil.hpp"


namespace OP2A{
namespace GRID{


class	FaceGeometry :public FaceGeoBasic, public FaceStencil
{
public:
	FaceGeometry()
	{

	}

	explicit FaceGeometry(const int ND, const FaceType f_type, bool extendedStencil) : FaceGeoBasic(ND, f_type), FaceStencil(extendedStencil)
	{

	}

	~FaceGeometry()
	{

	}

public:	// Grid Processing functions


};



}
}



#endif /* FACEGEOMETRY_HPP_ */
