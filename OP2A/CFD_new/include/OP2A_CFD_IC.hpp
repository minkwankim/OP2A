/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_IC.hpp
 * 			-  
 *  
 */

#ifndef OP2A_CFD_IC_HPP_
#define OP2A_CFD_IC_HPP_

#include <mkl.h>
#include <vector>
#include <limits>
#include "omp.h"

#include "../../DATA/include/OP2A_data.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"


using namespace std;

void Calculate_IC(CFD_variable_setup_ver2	variable_setup,
				vector< vector < double > > &rho_s,
				vector< vector < double > > &u,
				vector <double> &T, vector <double> &Tr, vector <double> &Tv, vector <double> &Te,
				vector< vector<double> > &Q, vector< vector<double> > &V, SPECIES_DATA_BASIC &species);
void Calculate_IC_ver2(CFD_variable_setup_ver2	variable_setup,	vector<double> &rho_s, vector<double> &u, vector<double> &T,  vector<double> &Q, vector<double> &V, SPECIES_DATA_BASIC &species);



void Calculate_IC_multi_fluid(vector<CFD_variable_setup_ver2>	&variable_setup,
								vector< vector < vector < double > > > &rho_s, vector< vector < vector < double > > > &u,
								vector< vector <double> > &T, vector< vector <double> > &Tr, vector< vector <double> > &Tv, vector< vector <double> > &Te,
								vector< vector< vector<double> > > &Q, vector< vector< vector<double> > > &V, vector<SPECIES_DATA_BASIC>	&species,
								int NUM_FLUID);

void Calculate_IC_multi_fluid_ver2(vector<SOL_CFD> &Solution,
									vector< vector < vector < double > > > &rho_s, vector< vector < vector < double > > > &u,
									vector< vector <double> > &T, vector< vector <double> > &Tr, vector< vector <double> > &Tv, vector< vector <double> > &Te,
									vector< vector< vector<double> > > &Q, vector< vector< vector<double> > > &V, vector<SPECIES_DATA_BASIC>	&species,
									int NUM_FLUID);


void Assign_IC (CFD_variable_setup_ver2	variable_setup, vector< vector<double> > &Q, vector< vector<double> > &V,
				SOL_CLASS_BASIC	&Cell_data_Q, SOL_CLASS_BASIC	&Cell_data_V, SPECIES_DATA_BASIC &species, int NCM, int ini_method, int NT);

void Assign_IC_ver2 (SOL_CFD &Solution, vector< vector<double> > &Q, vector< vector<double> > &V, int NCM, int ini_method, int NT);


void Initialize_flows_multi_fluids(vector<SOL_CFD> &Solution,
									vector< vector< vector<double> > > &Q,
									vector< vector< vector<double> > > &V,
									vector<SPECIES_DATA_BASIC>	&species,
									int NCM, int ini_method, int NUM_FLUID, int NT);

#endif /* OP2A_CFD_IC_HPP_ */
