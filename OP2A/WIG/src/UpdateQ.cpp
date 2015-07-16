/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 29, 2015
 *      			Author: Minkwan Kim
 *
 * UpdateQ.cpp
 * 			-  
 *  
 */



#include <omp.h>
#include "CFD/include/VariableSet.hpp"
#include "CFD/include/MiscFunctionsCFD.hpp"
#include "CFD/include/FluxInviscid.hpp"
#include "CFD/include/VariableSamplesWIG.hpp"
#include "Math/include/MathMisc.hpp"

#include "../include/OP2A_Application.hpp"


void ApplicationOP2A::UpdateQ()
{
	int indexE	= species_set.NS + grid.ND;

#pragma omp parallel for
	for (int c = 1; c <= grid.NCM; c++)
	{
#pragma ivdep
		for (int i = 0;  i<= grid.cells[c].data1D(indexQ).numData-1; i++)
		{
			grid.cells[c].data1D(indexQ)(i)	= grid.cells[c].data1D(indexQ)(i) + grid.cells[c].data1D(indexdQ)(i);
		}


		// Preliminary Error Check
		//	A. Density
#pragma ivdep
		for (int s = 0; s <= species_set.NS-1; s++)
		{
			/*
			if (grid.cells[c].data1D(indexQ)(s) < 0.0)
			{
				if (Math::fabs<double>(grid.cells[c].data1D(indexQ)(s)) <= 1.0e-16)
				{
					grid.cells[c].data1D(indexQ)(s) = 0.0;
				}
				else
				{
					std::ostringstream oss;
					oss << "It has negative species density: [Cell ID]: " << c << "  [Species ID]: " << s;
					throw Common::ExceptionNegativeValue (FromHere(), oss.str());
				}
			}
			*/

			if (Math::fabs<double>(grid.cells[c].data1D(indexQ)(s)) <= 1.0e-16)
			{
				grid.cells[c].data1D(indexQ)(s) = 0.0;
			}

			if (grid.cells[c].data1D(indexQ)(s) < 0.0)
			{
				std::ostringstream oss;
				oss << "It has negative species density: [Cell ID]: " << c << "  [Species ID]: " << s;
				throw Common::ExceptionNegativeValue (FromHere(), oss.str());
			}
		}

		// B. Energy
#pragma ivdep
		for (int m = indexE; m <= grid.cells[c].data1D(indexQ).numData-1; m++)
		{
			if (grid.cells[c].data1D(indexQ)(m) < 0.0)
			{
				std::ostringstream oss;
				oss << "It has negative energy: [Cell ID]: " << c << "  [Energy Mode]: " << m;
				throw Common::ExceptionNegativeValue (FromHere(), oss.str());
			}
		}
	}


#pragma omp parallel for
	for (int c = 1; c <= grid.NCM; c++)
	{
		CFD::VariableChange::From_Q(CFD_variabletype,
									grid.cells[c].data1D.data[indexQ],
									grid.cells[c].data1D.data[indexV],
									grid.cells[c].data1D.data[indexW],
									grid.cells[c].data1D.data[indexMIX],
									grid.cells[c].data1D.data[indexXs],
									grid.cells[c].data1D.data[indexYs],
									species_set, grid.ND, CFD_NT);
		// 1. Calculate V
		//CFD::VariableChange::Q_to_V(CFD_variabletype, CFD_NT, grid.cells[c].data1D(indexQ), species_set, grid.ND, grid.cells[c].data1D(indexV));

		// 2. Calculate W
		//CFD::VariableChange::V_to_W(CFD_variabletype, CFD_NT, grid.cells[c].data1D(indexV), species_set, grid.ND, grid.cells[c].data1D(indexW));
	}
	//Cell1DDataTreatement(1, true);
}
