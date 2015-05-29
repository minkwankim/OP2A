/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * transport_properties.hpp
 * 			-  
 *  
 */

#ifndef TRANSPORT_PROPERTIES_HPP_
#define TRANSPORT_PROPERTIES_HPP_


#include "../../CHEM/include/OP2A_chemistry.hpp"

double 	viscosity(double T, int model, Species* species_single);
double 	viscosity_collision(double T, double Te, double *gamms_s, Species *species_all, int NS);
double 	viscosity_mu_ver1(double T, int model, SPECIES &species);


double 	thermal_conductivity_Eucken(double mu_s, double T, Species* species_single, int MODE);
double 	diffusion_coeff_single_binary_Le_number(double rho, double k_tr, double Cp_tr, double Le);
double 	calculate_thermal_conductivity_Eucken(double mu_s, double T, SPECIES species, int MODE);
void 	calculate_thermal_conductivity(double mu_s, double Tv, SPECIES species, int method, vector<double> &kappa);




#define viscosity_mu	viscosity_mu_ver1

#endif /* TRANSPORT_PROPERTIES_HPP_ */
