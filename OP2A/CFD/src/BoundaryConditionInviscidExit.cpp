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


void BCInviscidImplicit::exitTypeBC(Data::DataStorage2D& J_plus, Data::DataStorage2D& J_minus, vector< vector<double> >& face_normal_vector, int NS, int ND, int NE)
{
	int VAR = NS + ND + NE;

#pragma ivdep
	for (int k = 0; k <= VAR-1; k++)
	{
		for (int l = 0; l <= VAR-1; l++)
		{
			J_plus(k,l)	+= J_minus(k,l);
			J_minus(k,l) = 0.0;
		}
	}
}



}
}
