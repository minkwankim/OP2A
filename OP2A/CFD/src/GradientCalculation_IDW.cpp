/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Sep 1, 2015
 *      			Author: Minkwan Kim
 *
 * GradientCalculation_IDW.cpp
 * 			-  
 *  
 */



#include <omp.h>
#include <time.h>

#include "CFD/include/GradientCalculation.hpp"

#include "Math/include/OP2A_Matrix.hpp"
#include "Math/include/MathMisc.hpp"
#include "Common/include/Time_Info.hpp"
#include "Common/include/MultiDimension.hpp"

namespace OP2A{
namespace CFD{


double GradientCalculation::GradientCalculation_IDW(vector<double>& x_f, vector< vector<double> >& x_i, vector <double>& phi_i, int N, int ND)
{
	vector< vector<double> > d_i	= Common::vector_2D<double>(N, ND, 0.0);

#pragma ivdep
	for (int i = 0; i <= N-1; i++)
	{
		for (int k = 0; k <= ND-1; k++)
		{
			d_i[i][k]	= x_i[i][k] - x_f[k];
		}
	}


	vector<double>	d_abs2(N, 0.0);
#pragma ivdep
	for (int i = 0; i <= N-1; i++)
	{
		for (int k = 0; k <= ND-1; k++)	d_abs2[i] += d_i[i][k]*d_i[i][k];
	}

	double aux_upper 	= 0.0;
	double aux_bottom	= 0.0;

	for (int i = 0; i <= N-1; i++)
	{
		aux_upper  += phi_i[i] / d_abs2[i];
		aux_bottom += 1.0 / d_abs2[i];
	}

	double phi_f	= aux_upper / aux_bottom;
	return (phi_f);
}


}
}
