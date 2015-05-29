/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: Dec 15, 2014
 *      			Author: Minkwan Kim
 *
 * global_variables.hpp
 * 			-  
 *  
 */

#ifndef GLOBAL_VARIABLES_HPP_
#define GLOBAL_VARIABLES_HPP_


extern int	*WHEREIS_CELL;				// Location of cells (cellptr/taskid)
extern int 	NUM_CELL_AT_TASK;			// Number of cells in the task
extern int	NP;							// Number of Processors
extern int	P;         					// Current Processor

extern int	OLD_OBJ;
extern int	NEW_OBJ;

extern bool	FLAG_SAMPLING;        		// Flag for sampling
extern int	FLAG_DIM;					// Flag for dimension (1 = Axisymmetric / 2 = 2D / 3 = 3D)
extern bool	FLAG_TERMINATION;      		// Flag for termination signal

extern double T_EPSIOLON;				// Relocation time tolerance



#ifndef MATH_ZERO
#define MATH_ZERO	1.0e-20
#endif

#endif /* GLOBAL_VARIABLES_HPP_ */
