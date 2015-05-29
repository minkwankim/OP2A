/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 15, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_reconstruct_limiter.cpp
 * 			-  
 *  
 */



#include "../include/OP2A_CFD_text.hpp"
#include "../include/OP2A_CFD_reconstruct.hpp"
#include "../../MATRIX/include/math_misc.hpp"

double minmod1(double r)
{
	double phi;

	phi = fmax(0.0, fmin (1.0, r));
	return phi;
}

double harmonic1(double r)
{
	double phi;

	phi = 4.0*r / pow(r+1.0, 2.0);
	return phi;
}

double superbee1(double r)
{
	double phi;
	phi	= fmax(0.0, fmin(2.0*r, 1.0));
	phi	= fmax(phi, fmin(r, 2.0));
	return phi;
}

double van_albada1(double r)
{
	double phi;
	phi = (r*r + r) / (r*r + 1.0);
	return phi;
}


double CFD_limiter(double r, double alpha, int method)
{
	double phi;
	switch (method)
	{
	case MINMOD:
		phi	= minmod1(r);
		break;

	case HARMONIC:
		phi	= harmonic1(r);
		break;

	case SUPERBEE:
		phi	= superbee1(r);
		break;

	case VAN_ALBADA:
		phi	= van_albada1(r);
		break;
	}

	if (phi <= 0.0)	phi	= 0.0;
	if (phi >= 2.0)	phi	= 2.0;

	phi	= phi / 2.0 * alpha;
	return phi;
}
