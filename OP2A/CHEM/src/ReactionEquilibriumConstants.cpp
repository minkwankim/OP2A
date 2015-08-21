/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 16, 2015
 *      			Author: Minkwan Kim
 *
 * ReactionEquilibriumConstants.cpp
 * 			-  
 *  
 */



#include <limits>

#include "Common/include/Exception_DimensionMatch.hpp"
#include "Common/include/Exception_NaNValue.hpp"
#include "Common/include/Exception_NegativeValue.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"
#include "Common/include/StringOps.hpp"

#include "Math/include/MathMisc.hpp"

#include "CHEM/include/ChemConstants.hpp"
#include "CHEM/include/Reaction.hpp"

using namespace std;

namespace OP2A{
namespace CHEM{





EquilibriumConstants::EquilibriumConstants():model(-1), m_num_data(0)
{

}

EquilibriumConstants::~EquilibriumConstants()
{

}

double EquilibriumConstants::Keq(double T, double n_mix)
{
	double a1	= Math::interpolate1(n, A1, n_mix);
	double a2	= Math::interpolate1(n, A2, n_mix);
	double a3	= Math::interpolate1(n, A3, n_mix);
	double a4	= Math::interpolate1(n, A4, n_mix);
	double a5	= Math::interpolate1(n, A5, n_mix);

	double log_keq;
	double Z	= 10000.0/T;

	switch (model)
	{
	case 85:
		log_keq	= a1 + a2*Z + a3*Z*Z + a4*pow(Z, 3.0) + a5*pow(Z, 4.0);
		break;

	case 90:
		log_keq	= a1/Z + a2 + a3*log(Z) + a4*Z + a5*Z*Z;
		break;

	case 94:
		log_keq	= a1 + a2*log(Z) + a3*Z + a4*pow(Z, 2.0) + a5*pow(Z, 3.0);
		break;
	}


	if (log_keq	!= log_keq)
	{
		throw Common::ExceptionNaNValue (FromHere(), "NaN value for Equilibrium constants");
	}

	if (log_keq == numeric_limits<double>::infinity())
	{
		throw Common::ExceptionInfiniteValue (FromHere(), "Infinite Value for Equilibrium constants");
	}

	return (exp(log_keq));
}



}
}
