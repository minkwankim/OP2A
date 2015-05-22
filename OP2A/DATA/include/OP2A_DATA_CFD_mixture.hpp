/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_DATA_CFD_mixture.hpp
 * 			-  
 *  
 */

#ifndef OP2A_DATA_CFD_MIXTURE_HPP_
#define OP2A_DATA_CFD_MIXTURE_HPP_

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

class CFD_mixture_data
{
public:
	double 	rho;
	double	Cv_bar;
	double 	Cv_bar_rot;
	double 	R_mix;
	double 	M_mix;
	double	gamma;

	vector<double> 	Xs;
	vector<double>	Ys;

	CFD_mixture_data();
	~CFD_mixture_data();

	void calculate_data(int NS, int NER, vector<double> Q, vector<double> Cv_tra_s, vector<double> Cv_rot_s, vector<double> R_s, vector<double> M_s);
};



#endif /* OP2A_DATA_CFD_MIXTURE_HPP_ */
