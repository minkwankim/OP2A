/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 16, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Time_integral.cpp
 * 			-  
 *  
 */




#include "../include/OP2A_CFD_Time_integration.hpp"


// Time integration
void CFD_time_integration(GRID_CLASS &grid, SOL_CFD &Solution, double dt, bool is_axis, unsigned int method, int nt)
{
	switch(method)
	{
	case EXPLICIT:
		CFD_time_integration_Explicit(grid, Solution, dt, is_axis, nt);
		break;

	case IMPLICIT_POINT:
		CFD_time_integration_Implicit_point(grid, Solution, dt, is_axis, nt);
		break;

	case IMPLICIT_LINE:
		//time_integration_Implicit_line(NVAR, grid, Qn, Rn, Jacobian_plus, Jacobian_minus, Jacobian_source, dQ, dt, is_axis, nt);
		break;
	}
}





