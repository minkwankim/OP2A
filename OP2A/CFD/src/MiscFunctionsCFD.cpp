/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 24, 2015
 *      			Author: Minkwan Kim
 *
 * MiscFunctions.CFD.cpp
 * 			-  
 *  
 */


#include "CFD/include/MiscFunctionsCFD.hpp"
#include "Math/include/MathMisc.hpp"



namespace OP2A{
namespace CFD{



	double CalculateCFLNumber(const int iter, const int iterBefore1, const double CFLStart, const double CFLMax)
	{
		double cfl;

		if (iter <= iterBefore1)	cfl	= Math::MK_fn_ver1(0, iterBefore1, CFLStart, 1.0, iter);
		else						cfl	= Math::MK_fn_ver1(iterBefore1, iterBefore1*10, 1.0, CFLMax, iter);

		return (cfl);
	}

	double Calculatedt(const double dX, const double U, const double Maxdt, const double CFL)
	{
		double dt = CFL * dX/U;

		if (dt <= 1.0e-14)	dt	= 1.0e-14;
		if (dt >= Maxdt)	dt	= Maxdt;

		return (dt);
	}



}
}
