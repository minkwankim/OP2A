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
	vector<double> relaxation(grid.NCM+1, 1.0);

#pragma omp parallel for
	for (int c = 1; c <= grid.NCM; c++)
	{
		vector<double> relaxation_factor(grid.cells[c].data1D(indexQ).numData, 1.0);
		Data::DataStorage Q_new("Q_new", grid.cells[c].data1D(indexQ).numData);
#pragma ivdep
		for (int i = 0;  i<= grid.cells[c].data1D(indexQ).numData-1; i++)
		{
			Q_new(i)	= grid.cells[c].data1D(indexQ)(i) + grid.cells[c].data1D(indexdQ)(i);
		}


		// Preliminary Error Check
		//	A. Density
#pragma ivdep
		for (int s = 0; s <= species_set.NS-1; s++)
		{
			if (Math::fabs<double>(Q_new(s)) <= 1.0e-16) Q_new(s) = 0.0;
			if (Q_new(s) < 0.0)
			{
				relaxation_factor[s] = - grid.cells[c].data1D(indexQ)(s) / grid.cells[c].data1D(indexdQ)(s);
			}
		}

		// B. Energy
#pragma ivdep
		for (int m = indexE; m <= grid.cells[c].data1D(indexQ).numData-1; m++)
		{
			if (Q_new(m) < 0.0)
			{
				relaxation_factor[m] = - grid.cells[c].data1D(indexQ)(m) / grid.cells[c].data1D(indexdQ)(m);
			}
		}

		relaxation[c] = Math::fmin<double>(relaxation_factor);
		if (relaxation[c] <= 0.0)
		{
			std::ostringstream oss;
			oss << "It has negative relaxation factor: " << relaxation[c] << "  [Cell ID]: " << c;
			throw Common::ExceptionNegativeValue (FromHere(), oss.str());
		}
	}

	double relaxation_min = Math::fmin<double>(relaxation);

	if (relaxation_min != 1.0)
	{
		relaxation_min *= 0.1;
		cout << "  [Relaxation Factor applied]: " << relaxation_min << endl;
	}

	int aaa;
	aaa = 0;

#pragma omp parallel for
	for (int c = 1; c <= grid.NCM; c++)
	{
#pragma ivdep
		for (int i = 0;  i<= grid.cells[c].data1D(indexQ).numData-1; i++)
		{
			grid.cells[c].data1D(indexQ)(i)	= grid.cells[c].data1D(indexQ)(i) + relaxation_min*grid.cells[c].data1D(indexdQ)(i);
		}

#pragma ivdep
		for (int s = 0; s <= species_set.NS-1; s++)
		{
			if (Math::fabs<double>(grid.cells[c].data1D(indexQ)(s)) <= 1.0e-16)
			{
				grid.cells[c].data1D(indexQ)(s) = 0.0;
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

/*

		// Preliminary Error Check
		//	A. Density
#pragma ivdep
		for (int s = 0; s <= species_set.NS-1; s++)
		{
			if (Math::fabs<double>(grid.cells[c].data1D(indexQ)(s)) <= 1.0e-16)
			{
				grid.cells[c].data1D(indexQ)(s) = 0.0;
			}

			if (grid.cells[c].data1D(indexQ)(s) < 0.0)
			{
				grid.cells[c].data1D(indexQ)(s) = 0.0;
				std::ostringstream oss;
				oss << "It has negative species density: [Cell ID]: " << c << "  [Species ID]: " << s;
				//throw Common::ExceptionNegativeValue (FromHere(), oss.str());
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

	*/
}
