/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 7, 2015
 *      			Author: Minkwan Kim
 *
 * cell_class.hpp
 * 			-  
 *  
 */

#ifndef CELL_CLASS_HPP_
#define CELL_CLASS_HPP_

#include <vector>
using namespace std;

#define CELL_CART_MAX_LEVEL	10

#include "../../DSMC/include/particle_object.hpp"













/*
 * Cell class for basic information
 */
class CELL_BASIC_2D
{
public:
	int	ID;
	unsigned int	type;					// type of cell
	unsigned int	BC;						// type of boundary condition

	unsigned int	NN;						// number of nodes in face
	vector<unsigned long int>	node;		// node IDs

	unsigned int	NF;						// number of faces in cell
	vector<unsigned long int>	face;		// face IDs


	double	S;								// Area
	double	characteristic_length;			// characteristic_length
	double 	dist_wall;						// distance to wall
	double 	x[2];							// position

	CELL_BASIC_2D();
	~CELL_BASIC_2D();

	// Internal functions
	void assign_NN_NF();
};

class CELL_BASIC_3D
{
public:
	unsigned int	type;					// type of cell
	unsigned int	BC;						// type of boundary condition

	unsigned int	NN;						// number of nodes in face
	vector<unsigned long int>	node;		// node IDs

	unsigned int	NF;						// number of faces in cell
	vector<unsigned long int>	face;		// face IDs


	double	S;								// Area
	double	characteristic_length;			// characteristic_length
	double 	dist_wall;						// distance to wall
	double 	x[3];							// position

	CELL_BASIC_3D();
	~CELL_BASIC_3D();

	// Internal functions
	void assign_NN_NF();
};






/*
 * Class for cell connectivity
 */
class CELL_CONNECTIVITY
{
public:
	unsigned int	N_neigh_C;				// number of neighboring cells
	vector<int>		neighbor;				// list of neighboring cells

	CELL_CONNECTIVITY();
	~CELL_CONNECTIVITY();
};


class CELL_CONNECTIVITY2
{
public:
	int	CN[2];
	int	CE[2];
	int	CW[2];
	int	CS[2];

	CELL_CONNECTIVITY2();
	~CELL_CONNECTIVITY2();
};





/*
 * Class for cartesian grid
 */
class CELL_CART_BASIC
{
public:
	unsigned int	cart_type;				// [0] Flow [1] Surface [2] Solid
	unsigned int	cut_in_cell_type;		// Cut-in-Cell type

	bool 			require_refinement;				// Flag for refinement
	bool			is_leaf;						// Flag for leaf

	unsigned int	mesh_level;				// mesh level
	bool			has_child;				// Refined

	CELL_CART_BASIC();
	~CELL_CART_BASIC();
};


class CELL_CART_INDEX_2D
{
public:
	unsigned int	IJK[2][CELL_CART_MAX_LEVEL];

	CELL_CART_INDEX_2D();
	~CELL_CART_INDEX_2D();
};

class CELL_CART_INDEX_3D
{
public:
	unsigned int	IJK[3][CELL_CART_MAX_LEVEL];

	CELL_CART_INDEX_3D();
	~CELL_CART_INDEX_3D();
};





#define CELL_BASIC		CELL_BASIC_2D
#define CELL_CART_INDEX	CELL_CART_INDEX_2D

#ifdef DIM_3D
#undef	CELL_BASIC
#define CELL_BASIC	CELL_BASIC_3D

#undef	CELL_CART_INDEX
#define CELL_CART_INDEX	CELL_CART_INDEX_3D
#endif







/*
 * BASIC class for CELL_CLASS
 */
class CELL_CLASS_BASIC	: public CELL_BASIC, public CELL_CONNECTIVITY
{
public:
	CELL_CLASS_BASIC();
	~CELL_CLASS_BASIC();
};


class CELL_CLASS_BASIC_CART	: public CELL_CLASS_BASIC, public CELL_CONNECTIVITY, public CELL_CART_BASIC, public CELL_CART_INDEX
{
public:
	CELL_CLASS_BASIC_CART();
	~CELL_CLASS_BASIC_CART();
};






/*
 * CELL CLASS for CFD simulation
 */
class CELL_CLASS_CFD_NORMAL	: public CELL_CLASS_BASIC
{
public:
	CELL_CLASS_CFD_NORMAL();
	~CELL_CLASS_CFD_NORMAL();
};

class CELL_CLASS_CFD_CART	: public CELL_CLASS_BASIC_CART
{
public:
	CELL_CLASS_CFD_CART				*parent;
	vector<CELL_CLASS_CFD_CART *>	child;

	CELL_CLASS_CFD_CART();
	~CELL_CLASS_CFD_CART();
};




/*
 * Cell class for DSMC/PIC simulation
 */
class CELL_CLASS_DSMC_NORMAL	: public CELL_CLASS_BASIC
{
public:
	unsigned long int	num_object[2];				// Number of objects (0:old / 1:new)
	particle_object		*object[2];					// Pointer to object in cell (0:old / 1:new)
	vector<int>			t_ratio;					// Ratio of time step



	CELL_CLASS_DSMC_NORMAL();
	~CELL_CLASS_DSMC_NORMAL();
};


class CELL_CLASS_DSMC_CART	: public CELL_CLASS_BASIC_CART
{
public:
	unsigned long int	num_object[2];				// Number of objects (0:old / 1:new)
	particle_object		*object[2];					// Pointer to object in cell (0:old / 1:new)
	vector<int>			t_ratio;					// Ratio of time step

	CELL_CLASS_DSMC_CART				*parent;
	vector<CELL_CLASS_DSMC_CART *>	child;

	CELL_CLASS_DSMC_CART();
	~CELL_CLASS_DSMC_CART();
};



#define	CELL_CLASS	CELL_CLASS_CFD_NORMAL

#ifdef	PARTICLE
#undef	CELL_CLASS
#define CELL_CLASS	CELL_CLASS_DSMC_NORMAL
#endif

#ifdef	CARTESIAN
#undef	CELL_CLASS
#define CELL_CLASS	CELL_CLASS_CFD_CART

#ifdef	PARTICLE
#undef	CELL_CLASS
#define CELL_CLASS	CELL_CLASS_DSMC_CART
#endif

#endif



#endif /* CELL_CLASS_HPP_ */
