/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 7, 2015
 *      			Author: Minkwan Kim
 *
 * face_class.hpp
 * 			-  
 *  
 */

#ifndef FACE_CLASS_HPP_
#define FACE_CLASS_HPP_

#include <vector>
#include "cell_class.hpp"
using namespace std;




/*
 * Basic information of face class
 */
class FACE_BASIC_3D
{
public:
	int				ID;						// ID of face
	unsigned int	type;					// type of face
	unsigned int	NN;						// number of nodes in face
	vector<unsigned long int>	node;		// node IDs
	double	n[3][3];						// [0] normal vector
											// [1] tangential vector
											// [2]tangential vector
	double S;								// Volume

	unsigned int	BC;						// Boundary condition index
	double 			x[3];					// position
	double 			dist_wall;				// distance to wall
	double 			n_dot_wall;

	FACE_BASIC_3D();
	~FACE_BASIC_3D();
};

class FACE_BASIC_2D
{
public:
	int				ID;						// ID of face
	unsigned int	type;					// type of face
	unsigned int	NN;						// number of nodes in face
	vector<unsigned long int>	node;		// node IDs
	double	n[2][2];						// [0] normal vector
											// [1] tangential vector
	double S;								// Length

	unsigned int	BC;						// Boundary condition index
	double 			x[2];					// position
	double 			dist_wall;				// distance to wall
	double 			n_dot_wall;

	FACE_BASIC_2D();
	~FACE_BASIC_2D();
};



/*
 * face class for stencill information
 */

class FACE_STENCIL
{
public:
	int cl[5];	// Left cells: 	[0] left		[1] left-left	[2] left-left-left		[3]cl1	[4]cl2
	int cr[5]; // Right cells:	[0]	right		[1] right-right	[2]	right-right-right	[3]cr1	[4]cr2


	FACE_STENCIL();
	~FACE_STENCIL();
};





/*
 * face class of additional information for CFD
 */
class FACE_ADD_CFD
{
public:
	double dist_to_wall;
	double n_dot_wall;

	FACE_ADD_CFD();
	~FACE_ADD_CFD();
};




#define FACE_BASIC	FACE_BASIC_2D

#ifdef DIM_3D
#undef	FACE_BASIC
#define FACE_BASIC	FACE_BASIC_3D
#endif



/*
 * FACE class for general
 */
class FACE_CLASS_NORMAL	: public FACE_BASIC, public FACE_STENCIL
{
public:
	FACE_CLASS_NORMAL();
	~FACE_CLASS_NORMAL();
};

class FACE_CLASS_CFD1	: public FACE_BASIC, public FACE_STENCIL
{
public:
	FACE_CLASS_CFD1();
	~FACE_CLASS_CFD1();
};

class FACE_CLASS_CFD2	: public FACE_BASIC, public FACE_STENCIL, public FACE_ADD_CFD
{
public:
	FACE_CLASS_CFD2();
	~FACE_CLASS_CFD2();
};

#define FACE_CLASS	FACE_CLASS_NORMAL

#ifdef CFD_NORMAL
#undef	FACE_CLASS
#define FACE_CLASS	FACE_CLASS_CFD1
#endif

#ifdef CFD_MSW
#undef	FACE_CLASS
#define FACE_CLASS	FACE_CLASS_CFD2
#endif

#endif /* FACE_CLASS_HPP_ */
