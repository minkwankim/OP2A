/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_DATA_CFD_transport_data.cpp
 * 			-  
 *  
 */

#ifndef OP2A_DATA_CFD_TRANSPORT_DATA_HPP_
#define OP2A_DATA_CFD_TRANSPORT_DATA_HPP_

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include "OP2A_DATA_CFD_setup.hpp"
#include "OP2A_DATA_CFD_transport_setup.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"

using namespace std;

class CFD_transport_properties_data
{
public:
	vector<double> 	diffusion_coefficient;

	vector<double> 	viscosity_coefficient;
	double 	viscosity_coefficient_mixture;

	vector<vector<double> >	thermal_conductivity;
	vector<double> 			thermal_conductivity_mixture;

	CFD_transport_properties_data();
	~CFD_transport_properties_data();

	void allocate_size(int ns, int ne);
	void mean_value(CFD_transport_properties_data &data1, CFD_transport_properties_data &data2, int NS, int NE);
	void calculate_transport_properties(CFD_variable_setup_ver2 &setup, vector<double> &V, vector<double> &Xs, SPECIES_DATA_BASIC &species, CFD_transport_properties_setup_ver2 viscosity_setup);
};



#endif /* OP2A_DATA_CFD_TRANSPORT_DATA_CPP_ */
