/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: Dec 12, 2014
 *      			Author: Minkwan Kim
 *
 * constants.hpp
 * 			-  
 *  
 */

#ifndef DSMC_CONSTANTS_HPP_
#define DSMC_CONSTANTS_HPP_

#define	MAX_NUM_SPECIES				5		//	Max. number of species
#define	MAX_NUM_BOUNDARY_SUMS		5		//	Max. number of boundary sums
#define MAX_NUM_OBJECT_PER_CELL		100000	//	Max. number of objects oer cell
#define	MAX_NUM_EDGES				30		//	Max. number of edges crossed per t_step

#define MATH_PI         			3.14159265358979323846
#define MATH_TWO_PI     			6.28318530717958647692
#define MATH_INFINITY   			3.40282347e+38F

#define CELL_ID_REMOVE				-100000


#define ERROR_MODULE				ERROR_MODULE_DSMC
#define ERROR_CODE_MOVE_PARTICLE	1


#endif /* CONSTANTS_HPP_ */
