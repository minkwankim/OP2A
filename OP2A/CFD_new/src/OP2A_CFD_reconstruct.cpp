/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 15, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_reconstruct.cpp
 * 			-  
 *  
 */



#include <cilk/cilk.h>
#include <omp.h>
#include <sstream>

#include "../include/OP2A_CFD_reconstruct.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../../SETUP/include/OP2A_setup_problem.hpp"
#include "../../MATRIX/include/OPPA_Matrix.hpp"

using namespace std;



/*
 * Reconstructor functions
 */
void reconstruct_MUSCL_ver3(vector<double>	&W_cll,	vector<double> &W_cl, vector<double> &W_cr, vector<double> &W_crr,
							vector<double>	&x_cll,	vector<double> &x_cl, vector<double> &x_cr, vector<double> &x_crr,
							vector<double>	&x_f, double kappa, int order, int limiter, int ND, int VAR,
							vector<double>	&W_L,	vector<double>	&W_R)
{


	// 1. Initialize Left and Right value
	W_L	= W_cl;
	W_R	= W_cr;



	// For high order
	if (order > 1)
	{
		double dx		= 0.0;
		double dx_m1	= 0.0;
		double dx_p1	= 0.0;
		double dx_f		= 0.0;
		double dx_fp	= 0.0;

		for (int k = 0; k <= ND-1; k++)
		{
			dx		+= pow((x_cr[k]		- x_cl[k]), 	2.0);
			dx_m1	+= pow((x_cl[k]		- x_cll[k]), 	2.0);
			dx_p1	+= pow((x_crr[k]	- x_cr[k]), 	2.0);

			dx_f	+= pow((x_f[k]		- x_cl[k]), 	2.0);
		}

		dx		= sqrt(dx);
		dx_m1	= sqrt(dx_m1);
		dx_p1	= sqrt(dx_p1);
		dx_f	= sqrt(dx_f);
		dx_fp	= dx - dx_f;


		// 3. W_L and W_R
		double 	r;
		double 	alpha;
		double	delta_W_m;
		double	delta_W;
		double	delta_W_p;

		double W_min, W_max;

		for (int k = 0; k <= VAR-1; k++)
		{
			delta_W	= (W_cr[k]	- W_cl[k]) 	/ dx;

			if (dx_m1 > 1.0e-10 && (W_cl[k]	- W_cll[k]) != 0.0)
			{
				delta_W_m	= (W_cl[k]	- W_cll[k]) / dx_m1;

				r			= delta_W / delta_W_m;
				alpha		= (dx/dx_f);

				W_L[k]	+=	CFD_limiter(r, alpha, limiter) * delta_W * dx_f;
			}
			else
			{
				W_L[k]	+=	delta_W * dx_f;
			}

			if (dx_p1 > 1.0e-10 && delta_W != 0.0)
			{
				delta_W_p	= (W_crr[k]	- W_cr[k]) / dx_p1;
				r			= delta_W_p / delta_W;
				alpha		= (dx / dx_fp);

				W_R[k]	-=	CFD_limiter(r, alpha, limiter) * delta_W * dx_fp;
			}
			else
			{
				W_R[k]	-= delta_W * dx_f;
			}

			W_min	= fmin(W_cl[k], W_cr[k]);
			W_max	= fmax(W_cl[k], W_cr[k]);

			if (W_L[k] < W_min)	W_L[k] = W_min;
			if (W_L[k] > W_max) W_L[k] = W_max;

			if (W_R[k] < W_min) W_R[k] = W_min;
			if (W_R[k] > W_max) W_R[k] = W_max;
		}
	}
}





