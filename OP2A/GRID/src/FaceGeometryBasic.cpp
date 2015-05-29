/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 29, 2015
 *      			Author: Minkwan Kim
 *
 * FaceGeometryBasic.cpp
 * 			-  
 *  
 */



#include "../include/FaceGeometryBasic.hpp"
#include "UTIL/include/OP2A_multi_dimension.hpp"

namespace OP2A{
namespace GRID{



FaceGeoBasic::FaceGeoBasic(const int ND, const FaceType f_type): type(f_type), x(ND, 0.0), n(ND, vector<double>(ND, 0.0))
{
	switch(type)
	{
	case FaceType::f_line:
		NN	= 2;
		break;

	case FaceType::f_quadrilateral:
		NN	= 4;
		break;

	case FaceType::f_triangle:
		NN = 3;
		break;
	}

	node_list.resize(NN);



	S 	= 0.0;
	BC	= BCType::interior;
	dist_wall	= 0.0;
	n_dot_wall	= 0.0;

	is_allocated = true;
}


FaceGeoBasic::FaceGeoBasic()
{
	type	= f_line;
	NN		= -1;
	S		= 0;

	BC		= interior;
	dist_wall 	= 0.0;
	n_dot_wall 	= 0.0;

	is_allocated = false;
}



FaceGeoBasic::~FaceGeoBasic()
{

}



void FaceGeoBasic::allocate(const int ND, const FaceType f_type)
{
	type	= f_type;

	x.resize(ND);
	n = vector_2D(ND, ND, 0.0);


	switch(type)
	{
	case FaceType::f_line:
		NN	= 2;
		break;

	case FaceType::f_quadrilateral:
		NN	= 4;
		break;

	case FaceType::f_triangle:
		NN = 3;
		break;
	}

	node_list.resize(NN);

	S 	= 0.0;
	BC	= BCType::interior;
	dist_wall	= 0.0;
	n_dot_wall	= 0.0;

	is_allocated = true;
}


bool FaceGeoBasic::allocation()
{
	return (is_allocated);
}





}
}
