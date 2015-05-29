/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 14, 2015
 *      			Author: Minkwan Kim
 *
 * OPPA_chemistry.hpp
 * 			-  
 *  
 */

#ifndef OPPA_CHEMISTRY_HPP_
#define OPPA_CHEMISTRY_HPP_


#include "Problem/include/OP2A_setup_constants.hpp"
#include "../../PHYSICS/include/physical_constants.hpp"
#include "OP2A_CHEM_constant.hpp"
#include "species.hpp"
#include "species_data_class.hpp"
#include "species_functions.hpp"
#include "reaction.hpp"


#define read_species		read_species_ver2
#define asg_recombination	asg_recombination_ver2

void chemical_production_rate(vector<SPECIES> &species_all, REACTION_DATA &reactions, vector<double> &rho_all, vector< vector<double> > &T_all, vector<double> &n_chem);
void calculate_reaction_rate(REACTION_DATA_ver2 &reaction_data, vector<double> &rhos, vector<double> &T, vector<double> &kf, vector<double> &kb, vector<double> &Rf, vector<double> &Rb);
void calculate_production_rate(REACTION_DATA_ver2 &reaction_data, vector<double> &Rf, vector<double> &Rb, vector<double> &S_chem);


#endif /* OPPA_CHEMISTRY_HPP_ */
