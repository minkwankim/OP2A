/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 9, 2015
 *      			Author: Minkwan Kim
 *
 * CalculateFluxInviscid.cpp
 * 			-  
 *  
 */



#include <omp.h>
#include "CFD/include/VariableSet.hpp"
#include "CFD/include/MiscFunctionsCFD.hpp"
#include "CFD/include/FluxInviscid.hpp"

#include "Math/include/MathMisc.hpp"

#include "../include/OP2A_Application.hpp"


void ApplicationOP2A::CalculateFluxInviscid()
{

	if (problem_setup.TIME_INTEGRATION_METHOD == 0)
	{
		CalculateFluxInviscidExplicit();
	}
	else
	{
		CalculateFluxInviscidImplicit();
		ApplyBCInviscidImplicit();
	}
}
