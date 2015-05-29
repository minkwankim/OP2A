/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 15, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_reconstruct.hpp
 * 			-  
 *  
 */

#ifndef OP2A_CFD_RECONSTRUCT_HPP_
#define OP2A_CFD_RECONSTRUCT_HPP_

#include <mkl.h>
#include <vector>
#include <limits>
#include "omp.h"


using namespace std;

double minmod1(double r);
double harmonic1(double r);
double superbee1(double r);
double van_albada1(double r);
double CFD_limiter(double r, double alpha, int method);


void reconstruct_MUSCL_ver3(vector<double>	&W_cll,	vector<double> &W_cl, vector<double> &W_cr, vector<double> &W_crr, vector<double>	&x_cll,	vector<double> &x_cl, vector<double> &x_cr, vector<double> &x_crr, vector<double>	&x_f, double kappa, int order, int limiter, int ND, int VAR, vector<double>	&W_L,	vector<double>	&W_R);




#endif /* OP2A_CFD_RECONSTRUCT_HPP_ */
