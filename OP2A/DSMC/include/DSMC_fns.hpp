/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: Dec 15, 2014
 *      			Author: Minkwan Kim
 *
 * DSMC_fns.hpp
 * 			-  
 *  
 */

#ifndef DSMC_FNS_HPP_
#define DSMC_FNS_HPP_

#include "constants.hpp"
#include "global_variables.hpp"
#include "particle_object.hpp"
#include "parallel.hpp"
#include "problem_setup.hpp"

#include "../../utilities/include/general_fns.hpp"
#include "../../mesh/include/cell_DSMC.hpp"

void move_particle(Cell_type **cell_p, particle_object *object, double &t_move, double &t_done, int &cell_ID_new, int &num_move, int &num_bad_move, int ND);
void calculate_cell(OPPA_DSMC_setup &problem, Cell_type	**cell_ptr, int ND, double &perf, double &perf2);





#endif /* DSMC_FNS_HPP_ */
