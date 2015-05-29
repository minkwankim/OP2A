/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 17, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_source_terms.hpp
 * 			-  
 *  
 */

#ifndef OP2A_CFD_SOURCE_TERMS_HPP_
#define OP2A_CFD_SOURCE_TERMS_HPP_

#include "../include/OP2A_CFD_Variable_change.hpp"
#include "../include/OP2A_CFD_reconstruct.hpp"
#include "../include/OP2A_CFD_Flux.hpp"


#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../../SETUP/include/OP2A_setup_problem.hpp"
#include "../../MATRIX/include/OPPA_Matrix.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../GRID/include/OP2A_grid.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"

void CFD_residue_source(GRID_CLASS &grid_data, SOL_CFD &Solutions, bool is_axis, int nt);


void CFD_source_chemical_NONEQ_ver2(vector <double *> &rho_s_ALL, vector<vector <double *> > &Ts_ALL, REACTION_DATA_ver2 &reactions, vector<double> &Rf, vector<double> &Rb, vector<double> &kf, vector<double> &kb, vector<double *> &S_chem, int nt);
void CFD_source_chemical_NONEQ_ver1(vector <double *> &rho_s_ALL, vector<vector <double *> > &Ts_ALL, REACTION_DATA_ver2 &reactions, vector<double *> &S_chem, int nt);
void CFD_source_chemical_NONEQ_ALL_explicit(int NCM, vector< vector <double *> > &rho_s_ALL, vector< vector<vector <double *> > > &Ts_ALL, REACTION_DATA_ver2 &reactions, vector< vector<double *> > &S_chem, int nt);
void CFD_source_chemical_NONEQ_ALL_implicit(int NCM, vector< vector <double *> > &rho_s_ALL, vector< vector<vector <double *> > > &Ts_ALL, REACTION_DATA_ver2 &reactions,  vector<vector<double> > &kf, vector<vector<double> > &Rf, vector<vector<double> > &kb, vector<vector<double> > &Rb, vector< vector<double *> > &S_chem, int nt);
void CFD_residue_CHEM_NONEQ(GRID_CLASS &grid_data, SOL_CFD &Solutions, bool is_axis, int nt);

void CFD_source_NONEQ_vib(GRID_CLASS &grid, SOL_CFD &Solution, SPECIES_DATA_BASIC &species, int nt);


#endif /* OP2A_CFD_SOURCE_TERMS_HPP_ */
