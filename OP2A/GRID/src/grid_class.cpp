/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 7, 2015
 *      			Author: Minkwan Kim
 *
 * grid_class.cpp
 * 			-  
 *  
 */

#include "../include/grid_class.hpp"

/*
 * BAsic information for grid class	: GRID_BASIC_INFO
 */
GRID_BASIC_INFO::GRID_BASIC_INFO()
{
	DIM	= 0;
	NNM	= 0;
	NFM	= 0;
	NCM	= 0;
	NGM = 0;
}

GRID_BASIC_INFO::~GRID_BASIC_INFO()
{

}







/*
 * Node data
 */
GRID_NODE::GRID_NODE()
{

}

GRID_NODE::~GRID_NODE()
{
	int nnm	= nodes.size();

	for (int n = 0; n <= nnm-1; n++)	nodes[n]	= NULL;
}

// Internal function
void GRID_NODE::resize(unsigned int nnm)
{
	whereis_node.resize(nnm);
	nodes.resize(nnm);
}






/*
 * FACE data
 */
GRID_FACE::GRID_FACE()
{

}

GRID_FACE::~GRID_FACE()
{
	int nfm	= faces.size();

	for (int n = 0; n <= nfm-1; n++)	faces[n]	= NULL;
}

// Internal function
void GRID_FACE::resize(unsigned int nfm)
{
	whereis_face.resize(nfm);
	faces.resize(nfm);
}



/*
 * CELL data
 */
GRID_CELL::GRID_CELL()
{

}

GRID_CELL::~GRID_CELL()
{
	int ncm	= cells.size();

	for (int n = 0; n <= ncm-1; n++)	cells[n]	= NULL;
}

// Internal function
void GRID_CELL::resize(unsigned int ncm)
{
	whereis_cell.resize(ncm);
	cells.resize(ncm);
}





/*
 * GRID_CLASS_DATA
 */
template <class    data_class>
GRID_CLASS_DATA<data_class>::GRID_CLASS_DATA()
{

}

template <class    data_class>
GRID_CLASS_DATA<data_class>::~GRID_CLASS_DATA()
{
	int num	= data_ptr.size();

	for (int n = 0; n <= num-1; n++)	data_ptr[n]	= NULL;
}

// Internal function
template <class    data_class>
void GRID_CLASS_DATA<data_class>::resize(unsigned int ncm)
{
	whereis.resize(ncm);
	//data_ptr.resize(ncm);
}







/*
 * ===========================================================================================
 * 				GRID DATA CLASS
 * 					- Version: 1.0
 * 					- DATE: Jan. 5. 2015
 * 					- BY Min Kwan Kim
 * ============================================================================================
 */
GRID_CLASS_ver1_5::GRID_CLASS_ver1_5()
{

}

GRID_CLASS_ver1_5::~GRID_CLASS_ver1_5()
{

}


// Internal Functions
void GRID_CLASS_ver1_5::resize(unsigned int nnm, unsigned int nfm, unsigned int ncm, unsigned int ncm_ghost)
{
	NNM	= ncm;
	NFM	= nfm;
	NCM	= ncm;
	NGM	= ncm_ghost;

	nodes.resize(NNM);
	faces.resize(NFM);
	cells.resize(NCM);
	cells_ghost.resize(NGM);
}


GRID_LINE::GRID_LINE()
{
	num_lines	= 0.0;
}

GRID_LINE::~GRID_LINE()
{

}


GRID_CLASS_ver1_6::GRID_CLASS_ver1_6()
{

}

GRID_CLASS_ver1_6::~GRID_CLASS_ver1_6()
{

}
