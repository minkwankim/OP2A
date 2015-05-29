/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: Dec 15, 2014
 *      			Author: Minkwan Kim
 *
 * problem_setup.cpp
 * 			-  
 *  
 */

#include <vector>
#include "../include/problem_setup.hpp"



OPPA_DSMC_setup::OPPA_DSMC_setup()
{
	intvl_save		= 100;			// Interval for writing restart file
	intvl_sample	= 100;			// Interval for taking sampling values
	intvl_eval		= 100;			// Interval for evaluating field values
	intvl_print		= 100;			// Interval for printing output
	intvl_domain	= 100;			// Interval for para domain decomposition

	n_current		= 0;			// Current simulation time step
	n_total			= 10000;        // Maximum number of time steps
	n_startup		= 100;        	// Number of time steps until stationary


	M_coll			= 0;			// Collision event counter
	M_bad_coll		= 0;    		// Bad collision event counter
	M_reac			= 0;			// Chemistry reaction event counter
	M_bad_reac		= 0;			// Bad Chemistry reaction event counters
	M_move			= 0;			// Movement event counters
	M_bad_move		= 0;    		// Bad Movement event counters

	W_particle		= 1.0;			// N_real / N_simulation
	t_ref			= 0.0;			// Reference simulation time step

	NCM				= 0;			// Number of cells in the task
	n_particle_tot	= 0;			// Total number of particles
}


OPPA_DSMC_setup::~OPPA_DSMC_setup()
{

}



GRID_INFO_DSMC::GRID_INFO_DSMC()
{
	NCM	= 0;
	ND	= 2;
}

GRID_INFO_DSMC::~GRID_INFO_DSMC()
{

}





/*
 * ====================
 * Internal functions
 * ====================
 */

// F- 01 reset_counter()
//			- Reset event counters
void OPPA_DSMC_setup::reset_counters()
{
	M_coll			= 0;			// Collision event counter
	M_bad_coll		= 0;    		// Bad collision event counter
	M_reac			= 0;			// Chemistry reaction event counter
	M_bad_reac		= 0;			// Bad Chemistry reaction event counters
	M_move			= 0;			// Movement event counters
	M_bad_move		= 0;    		// Bad Movement event counters
}


