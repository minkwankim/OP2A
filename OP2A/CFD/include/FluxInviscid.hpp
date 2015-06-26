/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 25, 2015
 *      			Author: Minkwan Kim
 *
 * FluxInviscid.hpp
 * 			-  
 *  
 */
#ifndef FLUXINVISCID_HPP_
#define FLUXINVISCID_HPP_



#include "CFD/include/VariableChange.hpp"
#include "CFD/include/VariableConstants.hpp"


namespace OP2A{
namespace CFD{

enum LimiterType
{
	MinMod		= 1,
	Harmonic 	= 2,
	VanAlbada	= 3,
	Superbee	= 4
};


class CFD_API Reconstruct: public Common::NonInstantiable< Reconstruct>
{
public:
	static double Limiter(double r, double alpha, int method);

	static void FirstOrder(Data::DataStorage& Wcl, Data::DataStorage& Wcr, Data::DataStorage& Wp, Data::DataStorage& Wm);
	static void SecondOrderMUSCL(Data::DataStorage& Wcll, Data::DataStorage& Wcl, Data::DataStorage& Wcr, Data::DataStorage& Wcrr,
								vector<double>&	xcll, vector<double>&	xcl, vector<double>&	xcr, vector<double>&	xcrr,
								vector<double>&	xf, bool use_limiter, int limiter,
								Data::DataStorage& Wp, Data::DataStorage& Wm);

};





class CFD_API FluxInviscid: public Common::NonInstantiable<FluxInviscid>
{
public:

};





}
}

#endif /* FLUXINVISCID_HPP_ */
