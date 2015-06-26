/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 23, 2015
 *      			Author: Minkwan Kim
 *
 * CalculateWallTemperature.cpp
 * 			-  
 *  
 */


#include "Math/include/MathMisc.hpp"
#include "Common/include/Exception_CatchInLoop.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"
#include "Common/include/Exception_NaNValue.hpp"

#include "CFD/include/BoundaryConditions.hpp"
#include "CFD/include/Variables.hpp"
#include "CHEM/include/ChemConstants.hpp"

namespace OP2A{
namespace CFD{

double CFD_CalculateWallTemperature_fn_wo_e(double Tw, double Tcl, double kappa_cl, double dn)
{
	double fnTw;

	fnTw	= STEFAN_BOLTZMANN_SI * pow(Tw, 4.0) - (kappa_cl/dn)*(Tcl - Tw);
	return (fnTw);
}

double CFD_CalculateWallTemperature_fn_w_e(double Tw, double Tcl, double kappa_cl, double dn, double Tref, vector<double> emissivity_below, vector<double> emissivity_above)
{
	double emissivity_ceoff;
	if (Tw <= Tref)	emissivity_ceoff	= emissivity_below[0] + emissivity_below[1]*Tw + emissivity_below[2]*pow(Tw, 2.0) + emissivity_below[3]*pow(Tw, 3.0) + emissivity_below[4]*pow(Tw, 4.0) + emissivity_below[5]*pow(Tw, 5.0);
	else			emissivity_ceoff	= emissivity_above[0] + emissivity_above[1]*Tw + emissivity_above[2]*pow(Tw, 2.0) + emissivity_above[3]*pow(Tw, 3.0) + emissivity_above[4]*pow(Tw, 4.0) + emissivity_above[5]*pow(Tw, 5.0);

	double fnTw;
	fnTw	= emissivity_ceoff * STEFAN_BOLTZMANN_SI * pow(Tw, 4.0) - (kappa_cl/dn)*(Tcl - Tw);
	return (fnTw);
}


double CFD_CalculateWallTemperature_dfn_wo_e(double Tw, double Tcl, double kappa_cl, double dn)
{
	double dfnTw;

	dfnTw	= 4.0 * STEFAN_BOLTZMANN_SI * pow(Tw, 3.0) + (kappa_cl/dn);
	return (dfnTw);
}

double CFD_CalculateWallTemperature_dfn_w_e(double Tw, double Tcl, double kappa_cl, double dn, double Tref, vector<double> emissivity_below, vector<double> emissivity_above)
{
	double emissivity_ceoff;
	if (Tw <= Tref)	emissivity_ceoff	= emissivity_below[0] + emissivity_below[1]*Tw + emissivity_below[2]*pow(Tw, 2.0) + emissivity_below[3]*pow(Tw, 3.0) + emissivity_below[4]*pow(Tw, 4.0) + emissivity_below[5]*pow(Tw, 5.0);
	else			emissivity_ceoff	= emissivity_above[0] + emissivity_above[1]*Tw + emissivity_above[2]*pow(Tw, 2.0) + emissivity_above[3]*pow(Tw, 3.0) + emissivity_above[4]*pow(Tw, 4.0) + emissivity_above[5]*pow(Tw, 5.0);

	double demissivity_ceoff;
	if (Tw <= Tref)	demissivity_ceoff	= emissivity_below[1] + 2.0*emissivity_below[2]*Tw + 3.0*emissivity_below[3]*pow(Tw, 2.0) + 4.0*emissivity_below[4]*pow(Tw, 3.0) + 5.0*emissivity_below[5]*pow(Tw, 4.0);
	else			demissivity_ceoff	= emissivity_above[1] + 2.0*emissivity_above[2]*Tw + 3.0*emissivity_above[3]*pow(Tw, 2.0) + 4.0*emissivity_above[4]*pow(Tw, 3.0) + 5.0*emissivity_above[5]*pow(Tw, 4.0);


	double dfnTw;
	dfnTw	= 4.0*emissivity_ceoff*STEFAN_BOLTZMANN_SI*pow(Tw, 3.0) + demissivity_ceoff*STEFAN_BOLTZMANN_SI*pow(Tw, 4.0) + (kappa_cl/dn);
	return (dfnTw);
}



double CFD_CalculateWallTemperature_wo_e(double Tcl, double kappa_cl, double dn, int iterMax, double eps0)
{
	double Tmax;
	double Tmin	= 100.0;

	double dT;
	double fn_new;
	double fn_old;
	double dfn_old;
	double Tnew;
	double Told;

	int iter;
	bool findValue;


	// 1. find Tmin /Tmax
	iter = 0;
	dT   = 100.0;
	Tnew = Tmin;
	Told = Tmin;
	findValue = false;

	fn_old = CFD_CalculateWallTemperature_fn_wo_e(Told, Tcl, kappa_cl, dn);
	while (findValue == false && iter <= 5000)
	{
		if (fn_old < 0.0)		Tnew += dT;
		else if (fn_old > 0.0)	Tnew -= dT;
		else					return (Told);

		if (Tnew <= 0.0)
		{
			Tnew	= 0.9*Tmin;
			dT		= 0.9*dT;
		}

		fn_new = CFD_CalculateWallTemperature_fn_wo_e(Tnew, Tcl, kappa_cl, dn);
		if (fn_old > 0)
		{
			if (fn_new > 0.0)
			{
				dT		= 1.1*dT;
				fn_old 	= fn_new;
			}
			else if (fn_new < 0.0)
			{
				Tmin = Tnew;
				Tmax = Told;
				findValue = true;
			}
			else
			{
				return (Tnew);
			}
		}
		else
		{
			if (fn_new > 0.0)
			{
				Tmin = Told;
				Tmax = Tnew;
				findValue = true;
			}
			else if (fn_new < 0.0)
			{
				dT		= 1.1*dT;
				fn_old 	= fn_new;
			}
			else
			{
				return (Tnew);
			}
		}

		iter++;
	}

	if (findValue == false)	throw Common::ExceptionCatchInLoop (FromHere(), "CAUGHT IN A LOOP: Twall Solver!");

	double signFn;

	fn_old = CFD_CalculateWallTemperature_fn_wo_e(Tmin, Tcl, kappa_cl, dn);
	signFn	= Math::fsgn<double>(fn_old);
	if (signFn > 0.0)throw Common::ExceptionCatchInLoop (FromHere(), "Problem in finding Tmin(Wall Temperature), fn(Tmin) should be smaller than 0");

	fn_old = CFD_CalculateWallTemperature_fn_wo_e(Tmax, Tcl, kappa_cl, dn);
	signFn	= Math::fsgn<double>(fn_old);
	if (signFn < 0.0)throw Common::ExceptionCatchInLoop (FromHere(), "Problem in finding Tmax(Wall Temperature), fn(Tmax) should be bigger than 0");



	// 2. Use Newton Method
	Told = Tmin;
	Tnew = Tmin;
	dT	 = 1.0;
	for (iter = 0; iter <= iterMax; iter++)
	{
		Told 	= Tnew;

		fn_old  = CFD_CalculateWallTemperature_fn_wo_e(Told, Tcl, kappa_cl, dn);
		dfn_old = CFD_CalculateWallTemperature_dfn_wo_e(Told, Tcl, kappa_cl, dn);
		dT 		= fn_old / dfn_old;

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

		if (Tnew > 1.1*Tmax) Tnew = Told + 10; 	//Tnew = Tmax - 2.0*Math::fabs<double>(dT);
		if (Tnew < 0.9*Tmin) Tnew = 0.9*Told; 	//Tmin + 0.1*Math::fabs<double>(dT);
	}
	cout << ".....FAILURE METHOD 1[NEWTON MEHTOD] TO Find Twall EQUATION. WILL TRY METHOD 2!!" << endl;



	// 3. Alternative Method
	dT 		= 10.0;
	Told 	= Tmin;

	fn_old  = CFD_CalculateWallTemperature_fn_wo_e(Told, Tcl, kappa_cl, dn);
	Tnew 	= Told + dT;

	for (iter = 0; iter <= 10*iterMax; iter++)
	{
		Told 	= Tnew;
		fn_new  = CFD_CalculateWallTemperature_fn_wo_e(Tnew, Tcl, kappa_cl, dn);

		signFn = Math::fsgn<double>(fn_new) * Math::fsgn<double>(fn_old);


		if (fn_old < 0.0)
		{
			if (fn_new > 0.0)
			{
				if (dT <= eps0)
				{
					return (Tnew);
				}
				else
				{
					dT = 0.7*dT;
					Tnew = Told - dT;
				}
			}
			else if (fn_new < 0.0)
			{
				dT = 1.2*dT;
				Tnew = Told + dT;
			}
			else
			{
				return (Tnew);
			}
		}
		else if (fn_old > 0)
		{
			if (fn_new > 0.0)
			{
				dT = 1.2*dT;
				Tnew = Told - dT;
			}
			else if (fn_new < 0.0)
			{
				if (dT <= eps0)
				{
					return (Tnew);
				}
				else
				{
					dT = 0.7*dT;
					Tnew = Told + dT;
				}
			}
			else
			{
				return (Tnew);
			}
		}
		else
		{
			return (Told);
		}


		Told = Tnew;
		if (Tnew > Tmax)	Tnew = Tmax - 0.1*Math::fabs<double>(dT);
		if (Tnew < Tmin)	Tnew = Tmin + 0.1*Math::fabs<double>(dT);
	}


	throw Common::ExceptionCatchInLoop (FromHere(), "Problem in finding Tv");
	return (0);
}

double CFD_CalculateWallTemperature_w_e(double Tcl, double kappa_cl, double dn, double Tref, vector<double> emissivity_below, vector<double> emissivity_above, int iterMax, double eps0)
{
	double Tmax;
	double Tmin	= 100.0;

	double dT;
	double fn_new;
	double fn_old;
	double dfn_old;
	double Tnew;
	double Told;

	int iter;
	bool findValue;


	// 1. find Tmin /Tmax
	iter = 0;
	dT   = 100.0;
	Tnew = Tmin;
	Told = Tmin;
	findValue = false;

	fn_old = CFD_CalculateWallTemperature_fn_w_e(Told, Tcl, kappa_cl, dn, Tref, emissivity_below, emissivity_above);

	while (findValue == false && iter <= 5000)
	{
		if (fn_old < 0.0)		Tnew += dT;
		else if (fn_old > 0.0)	Tnew -= dT;
		else					return (Told);

		if (Tnew <= 0.0)
		{
			Tnew	= 0.9*Tmin;
			dT		= 0.9*dT;
		}

		fn_new = CFD_CalculateWallTemperature_fn_w_e(Tnew, Tcl, kappa_cl, dn, Tref, emissivity_below, emissivity_above);
		if (fn_old > 0)
		{
			if (fn_new > 0.0)
			{
				dT		= 1.1*dT;
				fn_old 	= fn_new;
			}
			else if (fn_new < 0.0)
			{
				Tmin = Tnew;
				Tmax = Told;
				findValue = true;
			}
			else
			{
				return (Tnew);
			}
		}
		else
		{
			if (fn_new > 0.0)
			{
				Tmin = Told;
				Tmax = Tnew;
				findValue = true;
			}
			else if (fn_new < 0.0)
			{
				dT		= 1.1*dT;
				fn_old 	= fn_new;
			}
			else
			{
				return (Tnew);
			}
		}

		iter++;
	}

	if (findValue == false)	throw Common::ExceptionCatchInLoop (FromHere(), "CAUGHT IN A LOOP: Twall Solver!");

	double signFn;

	fn_old = CFD_CalculateWallTemperature_fn_w_e(Tmin, Tcl, kappa_cl, dn, Tref, emissivity_below, emissivity_above);
	signFn	= Math::fsgn<double>(fn_old);
	if (signFn > 0.0)throw Common::ExceptionCatchInLoop (FromHere(), "Problem in finding Tmin(Wall Temperature), fn(Tmin) should be smaller than 0");

	fn_old = CFD_CalculateWallTemperature_fn_w_e(Tmax, Tcl, kappa_cl, dn, Tref, emissivity_below, emissivity_above);
	signFn	= Math::fsgn<double>(fn_old);
	if (signFn < 0.0)throw Common::ExceptionCatchInLoop (FromHere(), "Problem in finding Tmax(Wall Temperature), fn(Tmax) should be bigger than 0");



	// 2. Use Newton Method
	Told = Tmin;
	Tnew = Tmin;
	dT	 = 1.0;
	for (iter = 0; iter <= iterMax; iter++)
	{
		Told 	= Tnew;

		fn_old  = CFD_CalculateWallTemperature_fn_w_e(Told, Tcl, kappa_cl, dn, Tref, emissivity_below, emissivity_above);
		dfn_old = CFD_CalculateWallTemperature_dfn_w_e(Told, Tcl, kappa_cl, dn, Tref, emissivity_below, emissivity_above);
		dT 		= fn_old / dfn_old;

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

		if (Tnew > 1.1*Tmax) Tnew = Told + 10; 	//Tnew = Tmax - 2.0*Math::fabs<double>(dT);
		if (Tnew < 0.9*Tmin) Tnew = 0.9*Told; 	//Tmin + 0.1*Math::fabs<double>(dT);
	}
	cout << ".....FAILURE METHOD 1[NEWTON MEHTOD] TO Find Twall EQUATION. WILL TRY METHOD 2!!" << endl;



	// 3. Alternative Method
	dT 		= 10.0;
	Told 	= Tmin;

	fn_old = CFD_CalculateWallTemperature_fn_w_e(Told, Tcl, kappa_cl, dn, Tref, emissivity_below, emissivity_above);
	Tnew 	= Told + dT;

	for (iter = 0; iter <= 10*iterMax; iter++)
	{
		Told 	= Tnew;

		fn_new = CFD_CalculateWallTemperature_fn_w_e(Tnew, Tcl, kappa_cl, dn, Tref, emissivity_below, emissivity_above);

		signFn = Math::fsgn<double>(fn_new) * Math::fsgn<double>(fn_old);


		if (fn_old < 0.0)
		{
			if (fn_new > 0.0)
			{
				if (dT <= eps0)
				{
					return (Tnew);
				}
				else
				{
					dT = 0.7*dT;
					Tnew = Told - dT;
				}
			}
			else if (fn_new < 0.0)
			{
				dT = 1.2*dT;
				Tnew = Told + dT;
			}
			else
			{
				return (Tnew);
			}
		}
		else if (fn_old > 0)
		{
			if (fn_new > 0.0)
			{
				dT = 1.2*dT;
				Tnew = Told - dT;
			}
			else if (fn_new < 0.0)
			{
				if (dT <= eps0)
				{
					return (Tnew);
				}
				else
				{
					dT = 0.7*dT;
					Tnew = Told + dT;
				}
			}
			else
			{
				return (Tnew);
			}
		}
		else
		{
			return (Told);
		}


		Told = Tnew;
		if (Tnew > Tmax)	Tnew = Tmax - 0.1*Math::fabs<double>(dT);
		if (Tnew < Tmin)	Tnew = Tmin + 0.1*Math::fabs<double>(dT);
	}


	throw Common::ExceptionCatchInLoop (FromHere(), "Problem in finding Tv");
	return (0);
}


double CFD_CalculateWallTemperature(double Tcl, double kappa_cl, double dn, bool use_emissivity, double Tref, vector<double> emissivity_below, vector<double> emissivity_above, int iterMax, double eps0)
{
	double Tw;

	if (use_emissivity == true)	Tw	= CFD_CalculateWallTemperature_w_e(Tcl, kappa_cl, dn, Tref, emissivity_below, emissivity_above, iterMax, eps0);
	else						Tw	= CFD_CalculateWallTemperature_wo_e(Tcl, kappa_cl, dn, iterMax, eps0);

	return (Tw);
}


}
}
