/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 29, 2015
 *      			Author: Minkwan Kim
 *
 * TimeIntegrate.cpp
 * 			-  
 *  
 */



#include <omp.h>
#include "CFD/include/VariableSet.hpp"
#include "CFD/include/MiscFunctionsCFD.hpp"
#include "CFD/include/FluxInviscid.hpp"

#include "Math/include/MathMisc.hpp"

#include "../include/OP2A_Application.hpp"

void ApplicationOP2A::TimeIntegrate()
{
	switch (problem_setup.TIME_INTEGRATION_METHOD)
	{
	case 0:
		TimeIntegrateExplicit();
		break;

	case 1:
		UpdateFluxJacobian();
		CalculateJacobianSourceTerm();
		/*
		CalculatedTdQ_dpdQ_at_cells();
		if (problem_setup.is_axisymmetric == true)	JacobianSourceAxis();
		*/

		TimeIntegrateImplicitPoint();
		break;
	}


	UpdateQ();
}
