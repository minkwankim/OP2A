/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 15, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_FLux.hpp
 * 			-  
 *  
 */

#ifndef OP2A_CFD_FLUX_HPP_
#define OP2A_CFD_FLUX_HPP_


#include "../include/OP2A_CFD_Variable_change.hpp"

#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../../SETUP/include/OP2A_setup_problem.hpp"
#include "../../MATRIX/include/OPPA_Matrix.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../GRID/include/OP2A_grid.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"


void CFD_Calculate_dTtra_dQ(CFD_variable_setup_ver2 setup, vector<double> &Q, vector<double> &V, CFD_mixture_data &variable_data, vector<double> &dT_dQ, SPECIES_DATA_BASIC &species);
void CFD_Calculate_dTrot_dQ(CFD_variable_setup_ver2 setup, vector<double> &Q, vector<double> &V, CFD_mixture_data &variable_data, vector<double> &dT_dQ, SPECIES_DATA_BASIC &species);
void CFD_Calculate_dTvib_dQ(CFD_variable_setup_ver2 setup, vector<double> &Q, vector<double> &V, CFD_mixture_data &variable_data, vector<double> &dT_dQ, SPECIES_DATA_BASIC &species);
void CFD_Calculate_dTe_dQ(CFD_variable_setup_ver2 setup, vector<double> &Q, vector<double> &V, CFD_mixture_data &variable_data, vector<double> &dT_dQ, SPECIES_DATA_BASIC &species);
void CFD_Calculate_dT_dQ(CFD_variable_setup_ver2 setup, vector<double> &Q, vector<double> &V, CFD_mixture_data &variable_data, vector< vector<double> > &dT_dQ, SPECIES_DATA_BASIC &species);

void CFD_Calculate_dp_dQ_type1(CFD_variable_setup_ver2 setup, vector<double> &V, CFD_mixture_data &variable_data, vector<vector<double> > &dT_dQ, vector<double> &dp_dQ,	SPECIES_DATA_BASIC &species);
void CFD_Calculate_dp_dQ_type2(CFD_variable_setup_ver2 setup, vector<double> &V, CFD_mixture_data &variable_data, vector<vector<double> > &dT_dQ, vector<double> &dp_dQ,	SPECIES_DATA_BASIC &species);
void CFD_Calculate_dp_dQ(CFD_variable_setup_ver2 setup, vector<double> &V, CFD_mixture_data &variable_data, vector<vector<double> > &dT_dQ, vector<double> &dp_dQ,	SPECIES_DATA_BASIC &species);
void CFD_Calculate_d2p_dQ2(CFD_variable_setup_ver2 setup, vector<double> &V, CFD_mixture_data &variable_data, vector < vector<double> > &d2p_dQ2,	SPECIES_DATA_BASIC &species);

void CFD_Calculate_da2_dQ(CFD_variable_setup_ver2 setup, vector<double> &Q, vector<double> &V, vector<double> &W, CFD_mixture_data &variable_data, vector<double> &dp_dQ, vector < vector<double> > &d2p_dQ2,	vector<double> &da2_dQ, SPECIES_DATA_BASIC &species);
void CFD_Calculate_da_dQ(CFD_variable_setup_ver2 setup, vector<double> &Q, vector<double> &V, vector<double> &W, CFD_mixture_data &variable_data, vector<double> &da2_dQ, double a, vector<double> &da_dQ, SPECIES_DATA_BASIC &species);
double CFD_Calculate_a2(CFD_variable_setup_ver2 setup, vector<double> &Q, vector<double> &W, CFD_mixture_data &variable_data, vector<double> &dp_dQ, SPECIES_DATA_BASIC &species);

void CFD_Flux_inviscid(GRID_CLASS &grid, SOL_CFD &Solution, SPECIES_DATA_BASIC &species, int order, int limiter, bool is_axis, double alpha, double x0, double eps0, int implicit, int nt);
void CFD_flux_1_2_FVM_single_ver3(CFD_variable_setup_ver2 setup, vector<double>	&Q_L, vector<double> &U_L, CFD_mixture_data &variable_data_L, vector<double>	&Q_R, vector<double> &U_R, CFD_mixture_data &variable_data_R, vector<double>	&F_1_2, vector< vector<double> > &n, double S, SPECIES_DATA_BASIC	&species, vector< vector<double> >	&Jacobian_plus, vector< vector<double> >	&Jacobian_minus, double dp, double dist_wall, double n_dot_wall, double alpha, double x0, double eps0);
void CFD_flux_1_2_FVM_single_ver2(CFD_variable_setup_ver2 setup, vector<double>	&Q_L, vector<double> &U_L, CFD_mixture_data &variable_data_L, vector<double>	&Q_R, vector<double> &U_R, CFD_mixture_data &variable_data_R, vector<double>	&F_1_2, vector< vector<double> > &n, double S, SPECIES_DATA_BASIC	&species, double dp, double dist_wall, double n_dot_wall, double alpha, double x0, double eps0);
void CFD_flux_1_2_AUSM_single(CFD_variable_setup_ver2 setup, vector<double>	&Q_L, vector<double> &U_L,  vector<double> &W_L, CFD_mixture_data &variable_data_L, vector<double>	&Q_R, vector<double> &U_R,  vector<double> &W_R, CFD_mixture_data &variable_data_R, vector<double>	&F_1_2, vector< vector<double> > &n, double S, SPECIES_DATA_BASIC	&species);


void CFD_Jacobian_inviscid_single(CFD_variable_setup_ver2 setup,	vector<double>	&Q, vector<double> &V, vector<double>	&W, CFD_mixture_data &variable_data, vector < vector<double> >	&dT_dQ, vector<double> &dp_dQ, vector < vector<double> >	&Jacob_inviscid, double signal, vector< vector<double> > &n, SPECIES_DATA_BASIC	&species, double dp, double dist_wall, double n_dot_wall, double x0, double eps0);


void CFD_residue_inviscid_ver2(GRID_CLASS &grid_data, SOL_CFD &Solution, bool is_axis, int nt);
void CFD_residue_inviscid_ver3(GRID_CLASS &grid_data, SOL_CFD &Solution, bool is_axis, int nt);



#endif /* OP2A_CFD_FLUX_HPP_ */
