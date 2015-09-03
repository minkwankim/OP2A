/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Aug 24, 2015
 *      			Author: Minkwan Kim
 *
 * species_collision_integral.cpp
 * 			-  
 *  
 */



#include <limits>

#include "Common/include/Exception_DimensionMatch.hpp"
#include "Common/include/Exception_NaNValue.hpp"
#include "Common/include/Exception_NegativeValue.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"

#include "CHEM/include/ChemConstants.hpp"
#include "CHEM/include/Species.hpp"
#include "Math/include/OP2A_Math.hpp"
#include "Math/include/MathMisc.hpp"



using namespace std;

namespace OP2A{
namespace CHEM{


double SpeciesTransport::piOmega1(const int l, const int r, const double T)
{
	double cross;
	cross = Math::interpolate1(Omega_temp[r][l], Omega[r][l], T);

	double aux;

	aux = MATH_PI * cross * 1.0e-20;

	return (aux);
}


double SpeciesTransport::piOmega0(const int l, const int r, const double T)
{
	double lnT	= log(T);

	double A = CCS[r][l][0];
	double B = CCS[r][l][1];
	double C = CCS[r][l][2];
	double D = CCS[r][l][3];

	double aux	= A* lnT*lnT + B*lnT + C;
	double collisionIntegral = D * pow(T, aux);
	collisionIntegral = collisionIntegral * 1.0e-20;	// Convert Angstron^2 to m^2

	return collisionIntegral;
}



double SpeciesTransport::piOmega2(const int l, const int r, const double T, const double ne, const double Te)	// Attractive force
{
	double Cn, cn, Dn;
	double Bs_a, Bs_b, Bs_c;
	double Cs_a, Cs_b, Cs_c;

	Bs_a = 1.0355;
	Bs_b = 0.2864;
	Bs_c = -4.4229;

	Cs_a = 0.3538;
	Cs_b = 0.1642;
	Cs_c = -3.8017;



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


	double debye_length_cgs;
	double Ts;
	double ne_cgs	= ne / 1.0e6;


	debye_length_cgs = sqrt(C_BOLTZMANN_CGS*Te / (4*MATH_PI*ne_cgs*ELECTRON_CHARGE_CGS_2));
	Ts				 = debye_length_cgs / (ELECTRON_CHARGE_CGS_2 / C_BOLTZMANN_CGS*T);

	double cross;
	cross = 5.0e15 * pow(debye_length_cgs/Ts, 2.0) * log(Dn*Ts * (1.0 - Cn*exp(-cn*Ts)) + 1);

	double collisionInt	= cross * MATH_PI * 1.0e-20;;

	return (collisionInt);
}



double SpeciesTransport::piOmega3(const int l, const int r, const double T, const double ne, const double Te)	// Repulsive force
{
	double Cn, cn, Dn;
	double Bs_a, Bs_b, Bs_c;
	double Cs_a, Cs_b, Cs_c;

	Bs_a = 1.0322;
	Bs_b = 0.4116;
	Bs_c = -3.9448;

	Cs_a = 0.3531;
	Cs_b = 0.2587;
	Cs_c = -3.1829;



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


	double debye_length_cgs;
	double Ts;
	double ne_cgs	= ne / 1.0e6;


	debye_length_cgs = sqrt(C_BOLTZMANN_CGS*Te / (4*MATH_PI*ne_cgs*ELECTRON_CHARGE_CGS_2));
	Ts				 = debye_length_cgs / (ELECTRON_CHARGE_CGS_2 / C_BOLTZMANN_CGS*T);

	double cross;
	cross = 5.0e15 * pow(debye_length_cgs/Ts, 2.0) * log(Dn*Ts * (1.0 - Cn*exp(-cn*Ts)) + 1);

	double collisionInt	= cross * MATH_PI * 1.0e-20;;

	return (collisionInt);
}



double SpeciesTransport::piOmega(const int l, const int r, const double T, const double ne, const double Te)
{
	double aux;
	switch (method_collision_integral[r])
	{
	case 0:
		aux = piOmega0(l, r, T);
		break;

	case 1:
		aux = piOmega1(l, r, T);
		break;

	case 2:
		aux = piOmega2(l, r, T, ne, Te);
		break;

	case 3:
		aux = piOmega3(l, r, T, ne, Te);
		break;
	}


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
