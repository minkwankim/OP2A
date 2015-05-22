/*
 *	Plasma Physics library - ver. 0.0
 *												by MINKWAN KIM
 *
 *				Last Modified Date: Sep 17, 2013
 *
 *	NOTE: IT CONTAINS FUNCTIONS IN THIS LIBRARY
 *
 *					Copyright (c) 2013 MINKWAN KIM
 *
 *  Created on: Jan 16, 2014
 *      Author: minkwan
 */

#ifndef PLASMA_PHYSICS_BASIC_PARAMETER_HPP_
#define PLASMA_PHYSICS_BASIC_PARAMETER_HPP_

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "./Plasma_constants.hpp"


// Debye length
double Debye_length_MKS(double ne, double Te);
double Debye_length_Gaussian(double ne, double Te);
double Debye_length_eV(double ne_mks, double Te_eV);

double V_thermal_MKS(double T, double m);
double Omega_p_MKS(double ne);
double Omega_p_CGS(double ne);



#endif /* PLASMA_PHYSICS_BASIC_PARAMETER_HPP_ */
