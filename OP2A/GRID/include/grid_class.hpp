/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 7, 2015
 *      			Author: Minkwan Kim
 *
 * grid_class.hpp
 * 			-  
 *  
 */

#ifndef GRID_CLASS_HPP_
#define GRID_CLASS_HPP_

#include <stddef.h>
#include <vector>

#include "node_class.hpp"
#include "face_class.hpp"
#include "cell_class.hpp"

using namespace std;


class GRID_BASIC_INFO
{
public:
	unsigned int	DIM;		// DIMENSION of grid
	unsigned int	NNM;		// number of nodes in grid
	unsigned int	NFM;		// number of faces in grid
	unsigned int	NCM;		// number of cells in grid
	unsigned int	NGM;		// number of ghost-cells in grid

	// Constructor and Destructor
	~GRID_BASIC_INFO();
	GRID_BASIC_INFO();
};




class GRID_NODE
{
public:
	vector<int>				whereis_node;	// location of nodes 		(node pointer/Task ID)
	vector<NODE_CLASS *>	nodes;

	GRID_NODE();
	~GRID_NODE();

	/*
	 * Internal functions
	 */
	void resize(unsigned int nnm);
};




class GRID_FACE
{
public:
	vector<int>				whereis_face;	// location of faces 		(face pointer/Task ID)
	vector<FACE_CLASS *>	faces;

	GRID_FACE();
	~GRID_FACE();

	/*
	 * Internal functions
	 */
	void resize(unsigned int nfm);
};




class GRID_CELL
{
public:
	vector<int>				whereis_cell;	// location of faces 		(face pointer/Task ID)
	vector<CELL_CLASS *>	cells;

	GRID_CELL();
	~GRID_CELL();

	/*
	 * Internal functions
	 */
	void resize(unsigned int ncm);
};




template <class    data_class>
class GRID_CLASS_DATA
{
public:
	vector<int>				whereis;	// location of data 		(data pointer/Task ID)
	vector<data_class *>	data_ptr;

	GRID_CLASS_DATA();
	~GRID_CLASS_DATA();

	/*
	 * Internal functions
	 */
	void resize(unsigned int ncm);
};






/*
 * ===========================================================================================
 * 				GRID DATA CLASS
 * 					- Version: 1.0
 * 					- DATE: Jan. 5. 2015
 * 					- BY Min Kwan Kim
 * ============================================================================================
 */
class GRID_CLASS_ver1	: public GRID_BASIC_INFO
{
public:
	GRID_NODE	NODES;
	GRID_FACE	FACES;
	GRID_CELL	CELLS;
	GRID_CELL	GHOSTS;
};


class GRID_CLASS_ver1_5	: public GRID_BASIC_INFO
{
public:
	GRID_CLASS_DATA <NODE_CLASS>	nodes;
	GRID_CLASS_DATA <FACE_CLASS>	faces;
	GRID_CLASS_DATA <CELL_CLASS>	cells;
	GRID_CLASS_DATA <CELL_CLASS>	cells_ghost;

	GRID_CLASS_ver1_5();
	~GRID_CLASS_ver1_5();

	// Internal Functions
	void	resize(unsigned int nnm, unsigned int nfm, unsigned int ncm, unsigned int ncm_ghost);
};



class GRID_LINE
{
public:
	int num_lines;
	vector <vector <int> >	lines;
	vector <vector <int> >	lines_bd;
	vector <int>			cell_line_info;


	GRID_LINE();
	~GRID_LINE();
};

class GRID_CLASS_ver1_6	: public GRID_CLASS_ver1_5
{
public:
	GRID_LINE	grid_line;

	GRID_CLASS_ver1_6();
	~GRID_CLASS_ver1_6();
};


#define GRID_CLASS	GRID_CLASS_ver1_6




#endif /* GRID_CLASS_HPP_ */
