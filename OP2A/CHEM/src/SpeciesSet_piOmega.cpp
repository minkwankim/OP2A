/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Aug 26, 2015
 *      			Author: Minkwan Kim
 *
 * SpeciesSet_piOmega.cpp
 * 			-  
 *  
 */

#include <limits>

#include "Common/include/Exception_DimensionMatch.hpp"
#include "Common/include/Exception_NaNValue.hpp"
#include "Common/include/Exception_NegativeValue.hpp"
#include "Common/include/Exception_OutOfRange.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"

#include "CHEM/include/ChemConstants.hpp"
#include "CHEM/include/Species.hpp"
#include "CHEM/include/SpeciesSet.hpp"

#include "Math/include/OP2A_Math.hpp"
#include "Math/include/MathMisc.hpp"


using namespace std;

namespace OP2A{
namespace CHEM{


// pi Omega: Method 0
double SpeciesSet::pi_Omega0(const int s, const int r, const int l, const double T)
{
	double lnT	= log(T);
	double aux	= CCS_A[s][r][l]* lnT*lnT + CCS_B[s][r][l]*lnT + CCS_C[s][r][l];
	double collisionIntegral = CCS_D[s][r][l] * pow(T, aux);

	collisionIntegral = collisionIntegral * 1.0e-20;	// Convert Angstron^2 to m^2

	return collisionIntegral;
}


// pi Omega: Method 1
double SpeciesSet::pi_Omega1(const int s, const int r, const int l, const double T)
{
	double cross;
	cross = Math::interpolate1(Omega_temp[s][r][l], Omega[s][r][l], T);

	double collisionIntegral = MATH_PI * cross * 1.0e-20;
	return (collisionIntegral);
}


// pi Omega: Method 2 (Attractive)
double SpeciesSet::pi_Omega2(const int s, const int r, const int l, const double T, const double ne, const double Te)
{
	double Cn, cn, Dn;

	switch (l)
	{
	case 1:
		Cn	= -0.476;
		cn	= 0.0313;
		Dn	= 0.784;
		break;

	case 2:
		Cn	= -0.146;
		cn	= 0.0377;
		Dn	= 1.262;
		break;
	}


	double ne_cgs 		= ne / 1.0e6;
	double Debye_length	= (C_BOLTZMANN_CGS*Te/(4.0*MATH_PI*ne_cgs*ELECTRON_CHARGE_CGS*ELECTRON_CHARGE_CGS));
	double Ts			= Debye_length / (ELECTRON_CHARGE_CGS*ELECTRON_CHARGE_CGS / (C_BOLTZMANN_CGS*T));

	double cross;
	cross = 5.0e15 * pow(Debye_length/Ts, 2.0) * log(Dn*Ts * (1.0 - Cn*exp(-cn*Ts)) + 1);

	double collisionIntegral	= cross * MATH_PI * 1.0e-20;;
	return (collisionIntegral);
}


// pi Omega: Method 3 (Repulsive)
double SpeciesSet::pi_Omega3(const int s, const int r, const int l, const double T, const double ne, const double Te)
{
	double Cn, cn, Dn;

	switch (l)
	{
	case 1:
		Cn	= 0.138;
		cn	= 0.0106;
		Dn	= 0.765;
		break;

	case 2:
		Cn	= 0.175;
		cn	= 0.0274;
		Dn	= 1.235;
		break;
	}


	double ne_cgs 		= ne / 1.0e6;
	double Debye_length	= (C_BOLTZMANN_CGS*Te/(4.0*MATH_PI*ne_cgs*ELECTRON_CHARGE_CGS*ELECTRON_CHARGE_CGS));
	double Ts			= Debye_length / (ELECTRON_CHARGE_CGS*ELECTRON_CHARGE_CGS / (C_BOLTZMANN_CGS*T));

	double cross;
	cross = 5.0e15 * pow(Debye_length/Ts, 2.0) * log(Dn*Ts * (1.0 - Cn*exp(-cn*Ts)) + 1);

	double collisionIntegral	= cross * MATH_PI * 1.0e-20;;
	return (collisionIntegral);
}



double SpeciesSet::pi_Omega(const int s, const int r, const int l, const double T, const double ne, const double Te)
{
	double aux;

	double Tp	= T;
	if (species[s].type == SpeciesType::Electron || species[r].type == SpeciesType::Electron) Tp = Te;


	switch (method_collision_integral[s][r])
	{
	case 0:
		aux = pi_Omega0(s, r, l, Tp);
		break;

	case 1:
		aux = pi_Omega1(s, r, l, Tp);
		break;

	case 2:
		aux = pi_Omega2(s, r, l, Tp, ne, Te);
		break;

	case 3:
		aux = pi_Omega3(s, r, l, Tp, ne, Te);
		break;
	}


	// CHECK ERROR
	if (aux	!= aux)
	{
		throw Common::ExceptionNaNValue (FromHere(), "NaN value for Collision integral ");
	}

	if (aux == numeric_limits<double>::infinity())
	{
		throw Common::ExceptionInfiniteValue (FromHere(), "Infinite Value for Collision integral");
	}

	if (aux < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Value for Collision integral");
	}

	return (aux);
}


double SpeciesSet::pi_Omega(const int s, const int r, const int l, const double T)
{
	double aux;

	if (method_collision_integral[s][r] == 0)
	{
		aux = pi_Omega0(s, r, l, T);
	}
	else if (method_collision_integral[s][r] == 0)
	{
		aux = pi_Omega1(s, r, l, T);
	}
	else
	{
		throw Common::ExceptionOutOfRange (FromHere(), "Collisional integral method: Need to provide ne and Te");
	}


	// CHECK ERROR
	if (aux	!= aux)
	{
		throw Common::ExceptionNaNValue (FromHere(), "NaN value for Collision integral ");
	}

	if (aux == numeric_limits<double>::infinity())
	{
		throw Common::ExceptionInfiniteValue (FromHere(), "Infinite Value for Collision integral");
	}

	if (aux < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Value for Collision integral");
	}

	return (aux);
}

}
}
