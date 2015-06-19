/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 19, 2015
 *      			Author: Minkwan Kim
 *
 * CalculateTve.cpp
 * 			-  
 *  
 */




#include <omp.h>

#include "CFD/include/CFDCommon.hpp"
#include "Math/include/MathMisc.hpp"
#include "Common/include/Exception_CatchInLoop.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"
#include "Common/include/Exception_NaNValue.hpp"



namespace OP2A{
namespace CFD{



double Fn_Tve(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, double Eve, double Tve, unsigned int CFD_NT)
{
	double rho_eve	= 0.0;
#pragma omp parallel for reduction(+:rho_eve)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		rho_eve	+= data_Q(s) * species_set.species[s].e_VE(Tve);
	}

	double o_result;
	o_result	= rho_eve - Eve;

	if (o_result == numeric_limits<double>::infinity())
	{
		throw Common::ExceptionInfiniteValue (FromHere(), "Inifnite Value for dFn_Tv");
	}

	if (o_result != o_result)
	{
		throw Common::ExceptionNaNValue (FromHere(), "NaN Value for dFn_Tv");
	}

	return(o_result);
}

double dFn_Tve(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, double Eve, double Tve, unsigned int CFD_NT)
{
	double rho_Cvve	= 0.0;
#pragma omp parallel for reduction(+:rho_Cvve)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		rho_Cvve	+= data_Q(s) * species_set.species[s].Cv_VE(Tve);
	}

	double o_result;
	o_result	= rho_Cvve;

	if (o_result == numeric_limits<double>::infinity())
	{
		throw Common::ExceptionInfiniteValue (FromHere(), "Inifnite Value for dFn_Tv");
	}

	if (o_result != o_result)
	{
		throw Common::ExceptionNaNValue (FromHere(), "Inifnite Value for dFn_Tv");
	}



	return(o_result);
}




double CFD_calculate_Tve(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, double Eve, unsigned int iterMax, double eps0, unsigned int CFD_NT)
{
	double Tmax;
	double Tmin	= 10.0;

	double dT;
	double fn_new;
	double fn_old;
	double dfn_old;
	double Tnew;
	double Told;

	int iter;
	bool findValue;

	// 1. find Tmax and Tmin
	iter = 0;
	dT   = 100.0;
	Tnew = Tmin;
	Told = Tmin;
	findValue = false;

	fn_old = Fn_Tve(data_Q, species_set, Eve, Told, CFD_NT);
	while (findValue == false && iter <= 5000)
	{
		if (fn_old < 0.0)	Tnew += dT;
		else				Tnew -= dT;

		if (Tnew <= 0.0)
		{
			Tnew = 0.9 * Tmin;
			dT = 0.9*dT;
		}

		fn_new = Fn_Tve(data_Q, species_set, Eve, Tnew, CFD_NT);

		if (fn_new*fn_old <= 0.0)
		{
			Tmax = Math::fmax<double> (Told, Tnew);
			Tmin = Math::fmin<double> (Told, Tnew);
			findValue = true;
		}
		else
		{
			dT = 1.1 * dT;
			fn_old = fn_new;
		}

		iter++;
	}
	if (findValue == false)	throw Common::ExceptionCatchInLoop (FromHere(), "CAUGHT IN A LOOP: Tve Solver!");


	double signFn;

	fn_old 	= Fn_Tve(data_Q, species_set, Eve, Tmin, CFD_NT);
	signFn	= Math::fsgn<double>(fn_old);
	if (signFn > 0.0)throw Common::ExceptionCatchInLoop (FromHere(), "Problem in finding Tmin, fn(Tmin) should be smaller than 0");

	fn_old 	= Fn_Tve(data_Q, species_set, Eve, Tmax, CFD_NT);
	signFn	= Math::fsgn<double>(fn_old);
	if (signFn < 0.0)throw Common::ExceptionCatchInLoop (FromHere(), "Problem in finding Tmax, fn(Tmax) should be bigger than 0");




	///////////////////////////////////////////////////////////////////////////////////////
	Told = Tmin;
	Tnew = Tmin;
	dT	 = 1.0;
	for (iter = 0; iter <= iterMax; iter++)
	{
		Told 	= Tnew;
		fn_old 	= Fn_Tve(data_Q, species_set, Eve, Told, CFD_NT);
		dfn_old = dFn_Tve(data_Q, species_set, Eve, Told, CFD_NT);

		dT = fn_old / dfn_old;

		if (dT != dT)									dT = -1.1*dT;
		if (dT == numeric_limits<double>::infinity())	dT = 100;
		if (dT >=  10000.0)								dT = 1000.0;
		if (dT <= -10000.0)								dT = -1000.0;


		if (Math::fabs<double>(dT) <= eps0)
		{
			return (Told);
		}
		else
		{
			Tnew	= Told - dT;
		}

		if (Tnew > 1.1*Tmax) Tnew = Told + 10; //Tnew = Tmax - 2.0*Math::fabs<double>(dT);
		if (Tnew < 0.9*Tmin) Tnew = 0.9*Told; //Tmin + 0.1*Math::fabs<double>(dT);
	}


	cout << ".....FAILURE METHOD 1[NEWTON MEHTOD] TO Find Tve EQUATION. WILL TRY METHOD 2!!" << endl;
	dT 		= 10.0;
	Told 	= Tmin;
	fn_old 	= Fn_Tve(data_Q, species_set, Eve, Told, CFD_NT);

	Tnew = Told + dT;
	for (iter = 0; iter <= 10*iterMax; iter++)
	{
		Told 	= Tnew;
		fn_new 	= Fn_Tve(data_Q, species_set, Eve, Tnew, CFD_NT);


		signFn = Math::fsgn<double>(fn_new) * Math::fsgn<double>(fn_old);

		if (signFn > 0.0)
		{
			dT = 1.2*dT;

			if (fn_new < 0.0)	Tnew = Told + dT;
			else				Tnew = Told - dT;
		}
		else if (signFn < 0.0)
		{
			dT = 0.7*dT;
			if (fn_new < 0.0)	Tnew = Told + dT;
			else				Tnew = Told - dT;
		}
		else
		{
			if (fn_new == 0.0)	return (Tnew);
		}

		Told = Tnew;

		if (Math::fabs<double>(Tnew - Told) <= eps0)	return (Tnew);


		if (Tnew > Tmax)	Tnew = Tmax - 0.1*Math::fabs<double>(dT);
		if (Tnew < Tmin)	Tnew = Tmin + 0.1*Math::fabs<double>(dT);
	}


	throw Common::ExceptionCatchInLoop (FromHere(), "Problem in finding Tv");
	return (0);
}




}
}
