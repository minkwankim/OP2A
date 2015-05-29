/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_problem_DSMC.hpp
 * 			-  
 *  
 */

#ifndef OP2A_PROBLEM_DSMC_HPP_
#define OP2A_PROBLEM_DSMC_HPP_



#include <string>
#include <vector>

using namespace std;


/*
 * ==============================
 * 		PROBLEM SETUP for DSMC
 * ==============================
 */
class PROBLEM_DSMC
{
public:
	unsigned int	itv_sample;				// Interval for taking sampling values
	unsigned int	itv_eval;				// Interval for evaluating field values
	unsigned int	itv_domain;				// Interval for para domain decomposition

	unsigned int	M_coll;					// Collision event counter
	unsigned int	M_bad_coll;    			// Bad collision event counter
	unsigned int 	M_reac;					// Chemistry reaction event counter
	unsigned int	M_bad_reac;				// Bad Chemistry reaction event counters
	unsigned int	M_move;					// Movement event counters
	unsigned int	M_bad_move;    			// Bad Movement event counters

	double			W_particle;			// N_real / N_simulation
	double			t_ref;				// Reference simulation time step

	unsigned int 	n_particle_tot;
	vector<double>	WEIGHT;		// Object weights for cells
	vector<double>	T_SCALE;	// Time scale factors for cells


	 PROBLEM_DSMC();
	~PROBLEM_DSMC();

	/*
	 * ====================
	 * Internal functions
	 * ====================
	 */

	// F- 01 reset_counter()
	//			- Reset event counters
	void reset_counters();
};


#endif /* OP2A_PROBLEM_DSMC_HPP_ */
