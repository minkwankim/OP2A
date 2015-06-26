/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 24, 2015
 *      			Author: Minkwan Kim
 *
 * CalculateCFL.cpp
 * 			-  
 *  
 */




#include "CFD/include/MiscFunctionsCFD.hpp"
#include "../include/OP2A_Application.hpp"



void ApplicationOP2A::CalcualteCFL()
{
	if (problem_setup.TIME_INTEGRATION_METHOD == 0)
	{
		CFLNumber	= problem_setup.CFL_start;
	}
	else
	{
		CFLNumber	= CFD::CalculateCFLNumber(iter, problem_setup.iteration_before_1, problem_setup.CFL_start, problem_setup.CFL_max);
	}
}
