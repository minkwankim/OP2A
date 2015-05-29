/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: Dec 15, 2014
 *      			Author: Minkwan Kim
 *
 * problem_setup.hpp
 * 			-  
 *  
 */

#ifndef PROBLEM_SETUP_HPP_
#define PROBLEM_SETUP_HPP_

#include <vector>
using namespace std;

/*
 * Data type for DSMC Grid
 */
class GRID_INFO_DSMC
{
public:
	int				NCM;		// Number of cells in the task
	int				ND;

	GRID_INFO_DSMC();
	~GRID_INFO_DSMC();
};


/*
 * Data type for DSMC Problem setup
 */
class OPPA_DSMC_setup
{
public:
	int	intvl_save;				// Interval for writing restart file
	int	intvl_sample;			// Interval for taking sampling values
	int	intvl_eval;				// Interval for evaluating field values
	int	intvl_print;			// Interval for printing output
	int	intvl_domain;			// Interval for para domain decomposition

	int	n_current;				// Current simulation time step
	int	n_total;            	// Maximum number of time steps
	int	n_startup;        		// Number of time steps until stationary



	int	M_coll;					// Collision event counter
	int	M_bad_coll;    			// Bad collision event counter
	int M_reac;					// Chemistry reaction event counter
	int	M_bad_reac;				// Bad Chemistry reaction event counters
	int	M_move;					// Movement event counters
	int	M_bad_move;    			// Bad Movement event counters

	double	W_particle;			// N_real / N_simulation
	double	t_ref;				// Reference simulation time step

	int				n_particle_tot;
	int				NCM;		// Number of cells in the task
	vector<double>	WEIGHT;		// Object weights for cells
	vector<double>	T_SCALE;	// Time scale factors for cells


	GRID_INFO_DSMC	grid;

	OPPA_DSMC_setup();
	~OPPA_DSMC_setup();


	/*
	 * ====================
	 * Internal functions
	 * ====================
	 */

	// F- 01 reset_counter()
	//			- Reset event counters
	void reset_counters();

};





#endif /* PROBLEM_SETUP_HPP_ */
