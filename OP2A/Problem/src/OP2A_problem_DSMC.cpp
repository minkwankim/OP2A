/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 13, 2015
 *      			Author: Minkwan Kim
 *
 * problem_DSMC.cpp
 * 			-  
 *  
 */


#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../include/OP2A_problem_DSMC.hpp"
#include "../include/OP2A_setup_constants.hpp"



 PROBLEM_DSMC::PROBLEM_DSMC()
{
	itv_sample	= 100;			// Interval for taking sampling values
	itv_eval	= 100;			// Interval for evaluating field values
	itv_domain	= 100;			// Interval for para domain decomposition

	M_coll			= 0;			// Collision event counter
	M_bad_coll		= 0;    		// Bad collision event counter
	M_reac			= 0;			// Chemistry reaction event counter
	M_bad_reac		= 0;			// Bad Chemistry reaction event counters
	M_move			= 0;			// Movement event counters
	M_bad_move		= 0;    		// Bad Movement event counters

	W_particle		= 1.0;			// N_real / N_simulation
	t_ref			= 0.0;			// Reference simulation time step

	n_particle_tot	= 0;			// Total number of particles
}


 PROBLEM_DSMC::~PROBLEM_DSMC()
{

}



/*
 * ====================
 * Internal functions
 * ====================
 */

// F- 01 reset_counter()
//			- Reset event counters
void PROBLEM_DSMC::reset_counters()
{
	M_coll			= 0;			// Collision event counter
	M_bad_coll		= 0;    		// Bad collision event counter
	M_reac			= 0;			// Chemistry reaction event counter
	M_bad_reac		= 0;			// Bad Chemistry reaction event counters
	M_move			= 0;			// Movement event counters
	M_bad_move		= 0;    		// Bad Movement event counters
}


