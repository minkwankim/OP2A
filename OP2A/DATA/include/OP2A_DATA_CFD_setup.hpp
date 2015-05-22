/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * CFD_variables.hpp
 * 			-  
 *  
 */

#ifndef DATA_CFD_VARIABLES_HPP_
#define DATA_CFD_VARIABLES_HPP_

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

class CFD_variable_setup_ver2
{
public:
	int	NS;
	int	ND;
	int NE;
	int VAR;

	int 				fluid_type;
	int 				energy_flag;
	vector <int> 		ID_T;

	CFD_variable_setup_ver2();
	~CFD_variable_setup_ver2();

	void assign_variables(int NER, int NEV, int NEE, int fluid_info);
};





#endif /* CFD_VARIABLES_HPP_ */
