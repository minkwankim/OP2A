/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 16, 2015
 *      			Author: Minkwan Kim
 *
 * ReactionRateCoeff.cpp
 * 			-  
 *  
 */

#include <limits>

#include "Common/include/Exception_DimensionMatch.hpp"
#include "Common/include/Exception_NaNValue.hpp"
#include "Common/include/Exception_NegativeValue.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"
#include "Common/include/StringOps.hpp"

#include "Math/include/OP2A_Math.hpp"

#include "CHEM/include/ChemConstants.hpp"
#include "CHEM/include/Reaction.hpp"

using namespace std;

namespace OP2A{
namespace CHEM{

Rate_coeff_Arrhenius::Rate_coeff_Arrhenius():Cf(0.0), nu(0.0), theta(0.0), k (-1.0)
{

}


Rate_coeff_Arrhenius::~Rate_coeff_Arrhenius()
{

}


double Rate_coeff_Arrhenius::ReactionRate(double Tc)
{
	double kf;

	kf = Cf * pow(Tc, nu) * exp(-theta *pow(Tc, k)) ;

	if (kf	!= kf)
	{
		throw Common::ExceptionNaNValue (FromHere(), "Nan value for Reaction coefficient ");
	}

	if (kf == numeric_limits<double>::infinity())
	{
		throw Common::ExceptionInfiniteValue (FromHere(), "Infinite Value for Reaction coefficient");
	}

	return (kf);
}


}
}
