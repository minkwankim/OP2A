/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 4, 2015
 *      			Author: Minkwan Kim
 *
 * GridConstant.hpp
 * 			-  
 *  
 */
#ifndef GRIDCONSTANT_HPP_
#define GRIDCONSTANT_HPP_



namespace OP2A{
namespace GRID{

#define	GRID_MAX_NUM_PATH		1000

enum GridDataType
{
	FLUENT	= 0
};

enum GridFluent
{
	FLU_INDEX_COMMENT  = 0,
	FLU_INDEX_DIMENSIONS  = 2,
	FLU_INDEX_NODE  = 10,
	FLU_INDEX_PERIODIC_SHADOW_FACES = 18,
	FLU_INDEX_CELL = 12,
	FLU_INDEX_FACE = 13,
	FLU_INDEX_BC = 45,
	FLU_INDEX_FACE_TREE = 59,
	FLU_INDEX_CELL_TREE = 58,
	FLU_INDEX_INTERFACE_FACE_PARENTS = 61,
};

enum BCType
{
	interior		= 0,
	wall			= 1,
	inlet			= 2,
	outlet			= 3,
	freestream		= 4,
	symmetric		= 5,
	axis			= 6,
	anode			= 7,
	cathode			= 8,
	dielectricwall	= 9
};

inline bool is_wall_typeBC(const int& bctype)
{
	if (bctype == BCType::wall || bctype == BCType::cathode || bctype == BCType::dielectricwall)	return (true);
	else return (false);
}


}
}


#endif /* GRIDCONSTANT_HPP_ */
