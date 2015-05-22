/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 8, 2015
 *      			Author: Minkwan Kim
 *
 * face_class.cpp
 * 			-  
 *  
 */



#include "../include/face_class.hpp"

/*
 * Basic information of face class
 */
FACE_BASIC_3D::FACE_BASIC_3D()
{
	ID	 = 0;
	type = 3;						// Type of face
	NN	 = 0;						// number of nodes in face
	node.reserve(4);				// node IDs
	S	= 0.0;

	BC	= 0;

	for (int i = 0; i <= 2; i++)	n[0][i] = 0.0;
	for (int i = 0; i <= 2; i++)	n[1][i] = 0.0;
	for (int i = 0; i <= 2; i++)	n[2][i] = 0.0;
	for (int i = 0; i <= 2; i++)	x[i] 	= 0.0;

	dist_wall	= 0.0;				// distance to wall
	n_dot_wall	= 0.0;
}

FACE_BASIC_3D::~FACE_BASIC_3D()
{

}


FACE_BASIC_2D::FACE_BASIC_2D()
{
	ID		= 0;
	type	= 2;					// type of face
	NN	= 0;						// number of nodes in face
	node.reserve(2);				// node IDs
	S	= 0.0;

	BC	= 0;

	for (int i = 0; i <= 1; i++)	n[0][i] = 0.0;
	for (int i = 0; i <= 1; i++)	n[1][i] = 0.0;
	for (int i = 0; i <= 1; i++)	x[i] 	= 0.0;

	dist_wall	= 0.0;				// distance to wall
	n_dot_wall	= 0.0;
}

FACE_BASIC_2D::~FACE_BASIC_2D()
{

}






/*
 * face class for stencill information
 */
FACE_STENCIL::FACE_STENCIL()
{
	for (int i = 0; i <= 4; i++)	cl[i]	= 0;
	for (int i = 0; i <= 4; i++)	cr[i]	= 0;
}

FACE_STENCIL::~FACE_STENCIL()
{

}




/*
* face class of additional information for CFD
*/
FACE_ADD_CFD::FACE_ADD_CFD()
{
	dist_to_wall	= 0.0;
	n_dot_wall		= 0.0;
}

FACE_ADD_CFD::~FACE_ADD_CFD()
{

}




/*
 * FACE class for general
 */
FACE_CLASS_NORMAL::FACE_CLASS_NORMAL()
{

}

FACE_CLASS_NORMAL::~FACE_CLASS_NORMAL()
{

}


FACE_CLASS_CFD1::FACE_CLASS_CFD1()
{

}

FACE_CLASS_CFD1::~FACE_CLASS_CFD1()
{

}


FACE_CLASS_CFD2::FACE_CLASS_CFD2()
{

}

FACE_CLASS_CFD2::~FACE_CLASS_CFD2()
{

}
