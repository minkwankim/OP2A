/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: May 29, 2014
 *      			Author: Minkwan Kim
 *
 * constant.hpp
 * 			-  
 *  
 */

#ifndef CONSTANT_GRID_HPP_
#define CONSTANT_GRID_HPP_





/*
 * =============================================
 * 	A. Grid Date reading
 * =============================================
 */
// 1. FLUENT grid
#define FLUENT			0
#define FLU_INDEX_COMMENT 0
#define FLU_INDEX_DIMENSIONS 2
#define FLU_INDEX_NODE 10
#define FLU_INDEX_PERIODIC_SHADOW_FACES 18
#define FLU_INDEX_CELL 12
#define FLU_INDEX_FACE 13
#define FLU_INDEX_BC 45
#define FLU_INDEX_FACE_TREE 59
#define FLU_INDEX_CELL_TREE 58
#define FLU_INDEX_INTERFACE_FACE_PARENTS 61

#define C_GHOST			-1
#define C_MIXED			0
#define C_TRIANGLE		1
#define C_TETRAHEDRON	2
#define C_QUADRILATERAL	3
#define C_HEXAHEDRON	4
#define C_PYRAMID		5
#define C_PRISM			6
#define C_WEDGE 		6


#define F_MIXED			0
#define F_LINE			2
#define F_TRIANGLE		3
#define F_QUADRILATERAL	4




/*
 * ====================================================
 * 	B. Boundary Condition Types
 * ====================================================
 */
#define MAX_ZONE_READ_MESH	1000

#define BC_INTERIOR			0
#define BC_WALL				1
#define BC_INLET			2
#define BC_OUTLET			3
#define BC_FREESTREAM		4
#define BC_SYMMETRY			5
#define BC_PERIODIC_SHADOW	6
#define BC_VELOCITY_INLET	7
#define BC_PERIODIC			8
#define BC_FAN				9
#define BC_MASS_FLOW_INLET	10
#define BC_INTERFACE		11
#define BC_PARENT			12
#define BC_OUTFLOW			13
#define BC_AXIS				14
#define BC_ANODE			15
#define BC_CATHODE			16
#define BC_DIELECTRIC_WALL	17
#define BC_LAST				BC_DIELECTRIC_WALL





#define N_MAX_SHAIRED_CELLS	8
#define F_MAX_NODES			4
#define C_MAX_NODES			8
#define C_MAX_FACES			6

#define	MAX_NUM_PATH		1000
#define	MAX_NUM_LINE		1000
#define	MAX_LINE_LENGTH		1000


#define		NORMAL		0
#define		TANGENTIAL1	1
#define		TANGENTIAL2	2

#define		OLD			0
#define		NEW			1

#define		CL			0
#define		CLL			1
#define		CLLL		2
#define		CL1			3
#define		CL2			4

#define		CR			0
#define		CRR			1
#define		CRRR		2
#define		CR1			3
#define		CR2			4

#define		MATH_ZERO_MESH	1.0e-20


#endif /* CONSTANT_HPP_ */
