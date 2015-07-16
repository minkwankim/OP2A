/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 10, 2015
 *      			Author: Minkwan Kim
 *
 * JacobianSourceTermAxisymmetric.cpp
 * 			-  
 *  
 */



#include <omp.h>
#include <time.h>

#include "CFD/include/Jacobians.hpp"
#include "Math/include/OP2A_Matrix.hpp"
#include "Math/include/MathMisc.hpp"
#include "Common/include/Time_Info.hpp"
#include "Common/include/MultiDimension.hpp"

namespace OP2A{
namespace CFD{

void Jacobians::Axisymmertic(Data::DataStorage& dpdQ, int NS, int ND, int NE, double S, Data::DataStorage2D& dSdQ)
{
	// 1. Initialize
	int VAR = NS + ND + NE;

#pragma ivdep
	for (int i = 0; i <= VAR-1; i++)
	{
		for (int j = 0; j <= VAR-1; j++)
		{
			dSdQ(i,j)	= 0.0;
		}
	}

	int indexP = NS + ND -1;
#pragma ivdep
	for (int i = 0; i <= VAR-1; i++)
	{
		dSdQ(indexP, i)	= dpdQ(i) * S;
	}
}




}
}
