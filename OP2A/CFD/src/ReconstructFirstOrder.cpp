/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 25, 2015
 *      			Author: Minkwan Kim
 *
 * ReconstructFirstOrder.cpp
 * 			-  
 *  
 */


#include "CFD/include/FluxInviscid.hpp"

namespace OP2A{
namespace CFD{


void Reconstruct::FirstOrder(Data::DataStorage& Wcl, Data::DataStorage& Wcr, Data::DataStorage& Wp, Data::DataStorage& Wm)
{
	Wp	= Wcl;
	Wm	= Wcr;
}








}
}
