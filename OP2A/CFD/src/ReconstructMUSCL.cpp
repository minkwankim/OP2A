/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 25, 2015
 *      			Author: Minkwan Kim
 *
 * ReconstructMUSCL.cpp
 * 			-  
 *  
 */

#include "CFD/include/FluxInviscid.hpp"
#include "Math/include/OP2A_Vector.hpp"
#include "Math/include/MathMisc.hpp"

#include "Common/include/Exception_InfiniteValue.hpp"
#include "Common/include/Exception_NaNValue.hpp"


namespace OP2A{
namespace CFD{

void Reconstruct:: SecondOrderMUSCL(Data::DataStorage& Wcll, Data::DataStorage& Wcl, Data::DataStorage& Wcr, Data::DataStorage& Wcrr,
								vector<double>&	xcll, vector<double>&	xcl, vector<double>&	xcr, vector<double>&	xcrr,
								vector<double>&	xf, bool use_limiter, int limiter,
								Data::DataStorage& Wp, Data::DataStorage& Wm)
{
	double dx, dx_plus, dx_minus, dx_f, dx_fp;

	Math::VECTOR	x(xcl, xcr);
	Math::VECTOR	x_plus(xcr, xcrr);
	Math::VECTOR	x_minus(xcll, xcl);
	Math::VECTOR	x_f(xcl, xf);

	dx	= x.length();
	dx_plus	= x_plus.length();
	dx_minus = x_minus.length();
	dx_f = x_f.length();
	dx_fp = dx - dx_f;

	if (dx_plus  == 0.0)	dx_plus		= 1.0E-14;
	if (dx_minus == 0.0)	dx_minus	= 1.0E-14;


	double dW;
	double dW_plus;
	double dW_minus;
	double r, alpha;
	double phi;

	double W_min;
	double W_max;


	for  (int index = 0; index <= Wcll.numData-1; index++)
	{
		dW			= (Wcr(index)	- Wcl(index))  / dx;
		dW_plus		= (Wcrr(index)	- Wcr(index))  / dx_plus;
		dW_minus	= (Wcl(index) 	- Wcll(index)) / dx_minus;

		if (dW_minus != 0.0 && use_limiter == true)
		{
			r		= dW / dW_minus;
			alpha	= dx / dx_f;
			phi		= Limiter(r, alpha, limiter);
		}
		else
		{
			phi = 1.0;
		}
		Wp(index)	= Wcl(index) + phi*dW*dx_f;
		if (Wp(index) != Wp(index))								throw Common::ExceptionNaNValue (FromHere(), "NaN value in ReconstructMUSCL");
		if (Wp(index) == numeric_limits<double>::infinity())	throw Common::ExceptionInfiniteValue (FromHere(), "Infinite value in ReconstructMUSCL");


		if (dW != 0.0 && use_limiter == true)
		{
			r		= dW_plus / dW;
			alpha	= dx / dx_fp;
			phi		= Limiter(r, alpha, limiter);
		}
		else
		{
			phi = 1.0;
		}

		Wm(index)	= Wcr(index) - phi*dW*dx_fp;
		if (Wm(index) != Wm(index))								throw Common::ExceptionNaNValue (FromHere(), "NaN value in ReconstructMUSCL");
		if (Wm(index) == numeric_limits<double>::infinity())	throw Common::ExceptionInfiniteValue (FromHere(), "Infinite value in ReconstructMUSCL");



		W_min	= Math::fmin<double>(Wcl(index), Wcr(index));
		W_max	= Math::fmax<double>(Wcl(index), Wcr(index));

		if (Wp(index) < W_min)	Wp(index)	= W_min;
		if (Wm(index) < W_min)	Wm(index)	= W_min;

		if (Wp(index) > W_max)	Wp(index)	= W_max;
		if (Wm(index) > W_max)	Wm(index)	= W_max;

	}
}



}
}
