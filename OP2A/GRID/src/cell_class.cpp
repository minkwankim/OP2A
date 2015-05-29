/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 8, 2015
 *      			Author: Minkwan Kim
 *
 * cell_class.cpp
 * 			-  
 *  
 */

#include <stddef.h>
#include "../include/cell_class.hpp"
#include "../include/constants_grid.hpp"

#include "../../UTIL/include/OP2A_utilities.hpp"

/*
 * Cell class for basic information
 */
CELL_BASIC_2D::CELL_BASIC_2D()
{
	ID	= 0;
	type	= 0;			// type of cell
	BC		= 0;			// type of boundary condition

	NN		= 0;			// number of nodes in face
	node.reserve(4);		// node IDs

	NF		= 0;			// number of faces in cell
	face.reserve(4);		// face IDs


	S						= 0.0;			// Area
	characteristic_length	= 0.0;			// characteristic_length
	dist_wall				= 0.0;			// distance to wall

	x[0]	= 0.0;
	x[1]	= 0.0;
}

CELL_BASIC_2D::~CELL_BASIC_2D()
{

}


// Internal functions
void CELL_BASIC_2D::assign_NN_NF()
{
	switch(type)
	{
	case C_TRIANGLE:
		NN	= 3;
		NF	= 3;
		break;

	case C_TETRAHEDRON:
		NN	= 4;
		NF	= 4;
		break;

	case C_QUADRILATERAL:
		NN	= 4;
		NF	= 4;
		break;

	case C_HEXAHEDRON:
		NN	= 8;
		NF	= 6;
		break;

	case C_PYRAMID:
		NN	= 5;
		NF	= 5;
		break;

	case C_WEDGE:
		NN	= 6;
		NF	= 5;
		break;

	default:
		Error_message_type	error;
		error.location_primary_name	= "cell_class.cpp";
		error.message = "It is not supported cell-type. Please use a different cell-type";
		error.print_message();
		break;
	}

	node.resize(NN);
	face.resize(NF);
}







CELL_BASIC_3D::CELL_BASIC_3D()
{
	type	= 0;			// type of cell
	BC		= 0;			// type of boundary condition

	NN		= 0;			// number of nodes in face
	node.reserve(8);		// node IDs

	NF		= 0;			// number of faces in cell
	face.reserve(6);		// face IDs


	S						= 0.0;			// Area
	characteristic_length	= 0.0;			// characteristic_length
	dist_wall				= 0.0;			// distance to wall

	x[0]	= 0.0;
	x[1]	= 0.0;
	x[2]	= 0.0;
}

CELL_BASIC_3D::~CELL_BASIC_3D()
{

}

void CELL_BASIC_3D::assign_NN_NF()
{
	switch(type)
	{
	case C_TRIANGLE:
		NN	= 3;
		NF	= 3;
		break;

	case C_TETRAHEDRON:
		NN	= 4;
		NF	= 4;
		break;

	case C_QUADRILATERAL:
		NN	= 4;
		NF	= 4;
		break;

	case C_HEXAHEDRON:
		NN	= 8;
		NF	= 6;
		break;

	case C_PYRAMID:
		NN	= 5;
		NF	= 5;
		break;

	case C_WEDGE:
		NN	= 6;
		NF	= 5;
		break;

	default:
		Error_message_type	error;
		error.location_primary_name	= "cell_class.cpp";
		error.message = "It is not supported cell-type. Please use a different cell-type";
		error.print_message();
		break;
	}

	node.resize(NN);
	face.resize(NF);
}




/*
 * Class for cell connectivity
 */
CELL_CONNECTIVITY::CELL_CONNECTIVITY()
{
	N_neigh_C	= 0;
	neighbor.reserve(8);
}

CELL_CONNECTIVITY::~CELL_CONNECTIVITY()
{

}


CELL_CONNECTIVITY2::CELL_CONNECTIVITY2()
{
	CN[0]	= 0;
	CN[1]	= 0;

	CE[0]	= 0;
	CE[1]	= 0;

	CW[0]	= 0;
	CW[1]	= 0;

	CS[0]	= 0;
	CS[1]	= 0;
}

CELL_CONNECTIVITY2::~CELL_CONNECTIVITY2()
{

}





/*
 * Class for cartesian grid
 */
CELL_CART_BASIC::CELL_CART_BASIC()
{
	cart_type			= 0;				// [0] Flow [1] Surface [2] Solid
	cut_in_cell_type	= 0;				// Cut-in-Cell type

	require_refinement	= false;			// Flag for refinement
	is_leaf				= true;				// Flag for leaf

	mesh_level			= 0;				// mesh level
	has_child			= false;				// Refined
}

CELL_CART_BASIC::~CELL_CART_BASIC()
{

}



CELL_CART_INDEX_2D::CELL_CART_INDEX_2D()
{
	for (int i = 0; i <= CELL_CART_MAX_LEVEL-1; i++)	IJK[0][i]	= 0;
	for (int i = 0; i <= CELL_CART_MAX_LEVEL-1; i++)	IJK[1][i]	= 0;
}

CELL_CART_INDEX_2D::~CELL_CART_INDEX_2D()
{

}


CELL_CART_INDEX_3D::CELL_CART_INDEX_3D()
{
	for (int i = 0; i <= CELL_CART_MAX_LEVEL-1; i++)	IJK[0][i]	= 0;
	for (int i = 0; i <= CELL_CART_MAX_LEVEL-1; i++)	IJK[1][i]	= 0;
	for (int i = 0; i <= CELL_CART_MAX_LEVEL-1; i++)	IJK[2][i]	= 0;
}

CELL_CART_INDEX_3D::~CELL_CART_INDEX_3D()
{

}




/*
 * BASIC class for CELL_CLASS
 */
CELL_CLASS_BASIC::CELL_CLASS_BASIC()
{

}

CELL_CLASS_BASIC::~CELL_CLASS_BASIC()
{

}


CELL_CLASS_BASIC_CART::CELL_CLASS_BASIC_CART()
{

}

CELL_CLASS_BASIC_CART::~CELL_CLASS_BASIC_CART()
{

}





/*
 * CELL CLASS for CFD simulation
 */
CELL_CLASS_CFD_NORMAL:: CELL_CLASS_CFD_NORMAL()
{

}

CELL_CLASS_CFD_NORMAL::~ CELL_CLASS_CFD_NORMAL()
{

}


CELL_CLASS_CFD_CART::CELL_CLASS_CFD_CART()
{
	parent	= NULL;
	child.reserve(8);
}

CELL_CLASS_CFD_CART::~CELL_CLASS_CFD_CART()
{
	parent	= NULL;
	for (int i = 0; i <= 7; i++)	child[i]	= NULL;
}





/*
 * Cell class for DSMC/PIC simulation
 */
CELL_CLASS_DSMC_NORMAL::CELL_CLASS_DSMC_NORMAL()
{
	num_object[0]	= 0;				// Number of objects (0:old / 1:new)
	num_object[1]	= 0;

	object[0]	= NULL;					// Pointer to object in cell (0:old / 1:new)
	object[1]	= NULL;
}

CELL_CLASS_DSMC_NORMAL::~CELL_CLASS_DSMC_NORMAL()
{
	object[0]	= NULL;					// Pointer to object in cell (0:old / 1:new)
	object[1]	= NULL;
}


CELL_CLASS_DSMC_CART::CELL_CLASS_DSMC_CART()
{
	num_object[0]	= 0;				// Number of objects (0:old / 1:new)
	num_object[1]	= 0;

	object[0]	= NULL;					// Pointer to object in cell (0:old / 1:new)
	object[1]	= NULL;

	parent	= NULL;
	child.reserve(8);
}

CELL_CLASS_DSMC_CART::~CELL_CLASS_DSMC_CART()
{
	object[0]	= NULL;					// Pointer to object in cell (0:old / 1:new)
	object[1]	= NULL;

	parent	= NULL;
	for (int i = 0; i <= 7; i++)	child[i]	= NULL;
}

