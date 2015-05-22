/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 15, 2015
 *      			Author: Minkwan Kim
 *
 * species_functions.hpp
 * 			-  
 *  
 */

#ifndef __SPECIES_FUNCTIONS_HPP_
#define __SPECIES_FUNCTIONS_HPP_


#include "OP2A_CHEM_constant.hpp"
#include "species.hpp"
#include "species_data_class.hpp"

#undef	SPECIES
#define SPECIES SPECIES_ver2

void read_species_ver2(vector<SPECIES> &species,	string file_name, int NS);
void asg_recombination_ver2(vector<SPECIES>& species, int NS);
void read_species_data_set(vector<SPECIES>& species, string file_name, int NS);







/*
 * Species Functions
 * 		ver - 1.0
 */
void FSpecies_read_species(Species *species, 			string file_name, int NS);
void FSpecies_read_species2(vector<Species>& species,	string file_name, int NS);


void FSpecies_asg_recombination(Species *species, int NS);
void FSpecies_asg_recombination2(vector<Species>& species, int NS);

void Read_species_DATA(Species *species, string file_name, int NS);
void Read_species_DATA2(vector<Species>& species, string file_name, int NS);
void Read_species_DATA_v3(vector<Species>& species, string file_name, int NS);


void FAsg_Omega_integral(Species* species, int NS);
void FAsg_Omega_integral2(vector<Species>& species, int NS);

void FSpecies_updata_Noneq_data(Species* species, int NS, int NONEQ_ROT, int NONEQ_VIB, int NONEQ_E);
void FSpecies_updata_Noneq_data2(vector<Species>& species, int NS, int NONEQ_ROT, int NONEQ_VIB, int NONEQ_E);


#endif /* SPECIES_FUNCTIONS_HPP_ */
