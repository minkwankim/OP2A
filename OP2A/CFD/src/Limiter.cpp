/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 25, 2015
 *      			Author: Minkwan Kim
 *
 * Limiter.cpp
 * 			-  
 *  
 */

#include <limits>
#include "CFD/include/FluxInviscid.hpp"
#include "Math/include/OP2A_Vector.hpp"
#include "Math/include/MathMisc.hpp"

#include "Common/include/Exception_OutOfRange.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"
#include "Common/include/Exception_NaNValue.hpp"



namespace OP2A{
namespace CFD{

double Reconstruct::Limiter(double r, double alpha, int method)
{
	double phi;

	if (r >= 0.0)
	{
		switch (static_cast<LimiterType>(method))
		{
		case LimiterType::MinMod:
			phi = Math::fmax<double>(0.0, Math::fmin<double>(1.0, r));
			break;

		case LimiterType::Harmonic:
			phi = 4.0*r / pow(r+1.0, 2.0);
			break;

		case LimiterType::Superbee:
			phi	= Math::fmax<double>(0.0, Math::fmin<double>(2.0*r, 1.0));
			phi	= Math::fmax<double>(phi, Math::fmin<double>(r, 2.0));
			break;

		case LimiterType::VanAlbada:
			phi = (r*r + r) / (r*r + 1.0);
			break;

		default:
			throw Common::ExceptionOutOfRange (FromHere(), "Out of range in CFD Limiter selection: The selected limited mode is not supported");
		}
	}
	else
	{
		phi = 0.0;
	}



	if (phi != phi)	throw Common::ExceptionNaNValue (FromHere(), "NaN value in CFD Limiter calculation");
	if (phi == numeric_limits<double>::infinity())	throw Common::ExceptionInfiniteValue (FromHere(), "Infinite value in CFD Limiter calculation");
	if (phi < 0.0)
	{
		throw Common::ExceptionOutOfRange (FromHere(), "Out of range in CFD Limiter calculation: Limiter should give positive value");
	}

	if (phi > 2.0)
	{
		throw Common::ExceptionOutOfRange (FromHere(), "Out of range in CFD Limiter calculation: Limited value should be less than 2");
	}


	phi	= phi / 2.0 * alpha;
	return phi;
}



}
}
