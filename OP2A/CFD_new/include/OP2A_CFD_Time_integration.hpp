/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 16, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Time_integration.hpp
 * 			-  
 *  
 */

#ifndef OP2A_CFD_TIME_INTEGRATION_HPP_
#define OP2A_CFD_TIME_INTEGRATION_HPP_

#include "../include/OP2A_CFD_Variable_change.hpp"

#include "../../UTIL/include/OP2A_utilities.hpp"
#include "Problem/include/OP2A_setup_problem.hpp"
#include "../../MATRIX/include/OPPA_Matrix.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../GRID/include/OP2A_grid.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"

void CFD_update_Q(GRID_CLASS &grid, SOL_CFD &Solution, int nt);

void CFD_time_integration(GRID_CLASS &grid, SOL_CFD &Solution, double dt, bool is_axis, unsigned int method, int nt);
void CFD_time_integration_Explicit(GRID_CLASS &grid, SOL_CFD &Solution, double dt, bool is_axis, int nt);
void CFD_time_integration_Implicit_point(GRID_CLASS &grid, SOL_CFD &Solution, double dt, bool is_axis, int nt);




#endif /* OP2A_CFD_TIME_INTEGRATION_HPP_ */
