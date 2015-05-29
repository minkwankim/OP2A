/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 16, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Jacobian.hpp
 * 			-  
 *  
 */

#ifndef OP2A_CFD_JACOBIAN_HPP_
#define OP2A_CFD_JACOBIAN_HPP_

#include "../include/OP2A_CFD_Variable_change.hpp"
#include "../include/OP2A_CFD_Flux.hpp"

#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../../SETUP/include/OP2A_setup_problem.hpp"
#include "../../MATRIX/include/OPPA_Matrix.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../GRID/include/OP2A_grid.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"



void CFD_Jacobian_dT_dQ_update_multi_fluid(int NCM, vector<CFD_variable_setup_ver2> &setup, vector<SOL_CFD> &Solutions, vector<SPECIES_DATA_BASIC> &species, REACTION_DATA_ver2 &reactions, vector< vector< vector<double> > > &dT_dQ, int num_fluid);

void CFD_Calculate_dRfb_dQ(double kfb, double Rfb, vector<double> &dkfb_dQ, vector<double> &alpha_j, vector<double> &rho_j, int NS, int VAR, vector<double> &dRfb_dQ);
void CFD_Calculate_dRfb_dQ_2(double kfb, double Rfb, vector<double> &dkfb_dQ, vector<double> &alpha_j, vector<double> &rho_j, int NS1, int ND1, int NE1, int VAR, vector<double> &dRfb_dQ);
void CFD_Calculate_dkf_dQ(double kf, double Tc, double nu_k, double theta_k, vector<double> &dTc_dQ, int NVAR, vector<double>	&dkf_dQ);
void CFD_Calculate_dkb_dQ(double kb, double Tbc, double nu_k, double theta_k, double dKeq_dT_over_Keq, vector<double> &dTbc_dQ, int NVAR, vector<double>	&dkb_dQ);
void CFD_Jacobian_chemistry_single_fluid(GRID_CLASS &grid, SOL_CFD &Solution, vector< vector<double> > &kf, vector< vector<double> > &Rf, vector< vector<double> > &kb, vector< vector<double> > &Rb, REACTION_DATA_ver2 &reactions, vector< vector< vector<double> > > &dS_chem_dQ, bool is_axis, int nt);



void CFD_Jacobian_axis_inviscid(CFD_variable_setup_ver2 setup,	vector<double> &dp_dQ, vector< vector<double> >	&Jacob_axis_inviscid, double S);
void CFD_Jacobian_Source_terms(GRID_CLASS &grid, SOL_CFD &Solution,	SPECIES_DATA_BASIC &species, REACTION_DATA_ver2 &reactions, vector<vector<double> > &kf, vector<vector<double> > &Rf, vector<vector<double> > &kb, vector<vector<double> > &Rb, bool is_axis, bool is_viscous, bool is_chemical, int nt);


#endif /* OP2A_CFD_JACOBIAN_HPP_ */
