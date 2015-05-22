/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * mixture_properties.hpp
 * 			-  
 *  
 */

#ifndef OP2A_MIXTURE_PROPERTIES_HPP_
#define OP2A_MIXTURE_PROPERTIES_HPP_



#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector>


using namespace std;

/*
 * =====================================
 * 		Mixture properties functions
 * ======================================
 */
double 	cal_mixture_properties_ver1(vector<double>Q, int NS, vector<double> Yi, vector<double> Xi, int mode);
double 	cal_mixing_Wilke_main_algorithm(vector<double> Qs, vector<double> Xs, vector<double> phi_s, unsigned int NS);
void 	mixing_Wilke_phi_ver1(vector<double> mu_s, vector<double> Xs, vector<double> Ms, vector<double> &phi_s, unsigned int NS);
double 	mixing_Wilke(vector<double> Qs, vector<double> mu_s, vector<double> Xs, vector<double> Ms, unsigned int NS);


#define cal_mixture_properties		cal_mixture_properties_ver1
#define mixing_Wilke_phi		mixing_Wilke_phi_ver1

#endif /* MIXTURE_PROPERTIES_HPP_ */
