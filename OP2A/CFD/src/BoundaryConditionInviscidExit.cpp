/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 23, 2015
 *      			Author: Minkwan Kim
 *
 * BoundaryConditionInviscidExit.hpp
 * 			-  
 *  
 */


#include "CFD/include/BoundaryConditions.hpp"

namespace OP2A{
namespace CFD{

void BCInviscid::exitTypeBC(Data::DataStorage& Qcl, Data::DataStorage& Qcr, int ND)
{
	int numVar	= Qcl.numData;

//#pragma omp parallel for
#pragma ivdep
	for (int i = 0; i <= numVar-1; i++)	Qcr(i)	= Qcl(i);
}



}
}
