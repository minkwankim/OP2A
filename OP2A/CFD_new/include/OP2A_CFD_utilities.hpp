/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_utilities.hpp
 * 			-  
 *  
 */

#ifndef OP2A_CFD_UTILITIES_HPP_
#define OP2A_CFD_UTILITIES_HPP_


#include <mkl.h>
#include <vector>
#include <limits>
#include "omp.h"

#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../../GRID/include/OP2A_grid.hpp"
#include "../../DATA/include/OP2A_data.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"


double calculate_dt(int NCM, vector<double> &dX, vector<double> &U, double dt_MAX, double CFL);
double calculate_dt_multi(int NCM, vector<double> &dX, vector< vector<double> > &U, double max_dt, double CFL, int NUM_FLUID);



int	energy_mode_flag(int NER, int NEV, int NEE);
void temperature_mode_table(int flag, vector<int> &ID_T);


void OP2A_CFD_assign_rho_and_T_for_all_fluid(int NNM, vector<SOL_CFD> &Solutions, vector<SPECIES_DATA_BASIC>	&species, vector<vector <double *> > &rho_s_ALL, vector<vector<vector <double *> > > &Ts_ALL, int num_fluid);
void OP2A_CFD_assign_rho_and_T_for_all_fluid_ver2(int NCM, vector<SOL_CFD> &Solutions, vector<SPECIES_DATA_BASIC>	&species, vector<vector <double> > &rho_s_ALL, vector<vector<vector <double> > > &Ts_ALL, int num_fluid);

void calculate_node_value_find_make_lists(GRID_CLASS &grid, vector < vector<int> >	&node_shared_cell_list, vector < vector<double> >	&node_shared_cell_weighting, int nt);
void calculate_node_value(GRID_CLASS &grid, SOL_CLASS_DATA &cell_Q,  SOL_CLASS_DATA &cellg_Q,  SOL_CLASS_DATA &node_Q, vector < vector<int> >	&node_shared_cell_list, vector < vector<double> >	&node_shared_cell_weighting, vector<int> variable_index, int num_variable, int nt);

void CFD_calc_residue_norms(int NS, int ND, int NE, GRID_CLASS &grid_data, vector< vector<double> > &Residue_c, double &RHS_max, int &RHS_max_i, double &RHS2, int nt);


#endif /* OP2A_CFD_UTILITIES_HPP_ */
