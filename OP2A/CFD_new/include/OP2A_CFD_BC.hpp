/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 15, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_BC.hpp
 * 			-  
 *  
 */

#ifndef OP2A_CFD_BC_HPP_
#define OP2A_CFD_BC_HPP_



#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include "../include/OP2A_CFD_Variable_change.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../../GRID/include/OP2A_grid.hpp"


void CFD_BC_Inviscid_wall(CFD_variable_setup_ver2 setup, CELL_CLASS &cell_cl, FACE_CLASS &face, CELL_CLASS &cell_cr, vector<double> &Q_cl, vector<double> &Q_cr);
void CFD_BC_Inviscid_inlet(CFD_variable_setup_ver2 setup, CELL_CLASS &cell_cl, CELL_CLASS &cell_cr, vector<double> &Q_cl, vector<double> &Q_cr, vector<double> &Q_inlet);
void CFD_BC_Inviscid_exit(CFD_variable_setup_ver2 setup, CELL_CLASS &cell_cl, CELL_CLASS &cell_cr, vector<double> &Q_cl, vector<double> &Q_cr);
void CFD_assign_BC_inviscid(CFD_variable_setup_ver2 setup, vector < vector <double> > &Q_inlet, GRID_CLASS &grid, vector < vector<double> > &Qc, vector < vector<double> > &Qgc, int nt);
void CFD_assign_BC_inviscid_complete(SOL_CFD &Solution_data, GRID_CLASS &grid, vector < vector <double> > &Q_IC, vector < vector <double> > &V_IC, SPECIES_DATA_BASIC &species, int nt);
void CFD_assign_BC_inviscid_implicit(CFD_variable_setup_ver2 setup, GRID_CLASS &grid, vector< vector< vector<double> > > &Jacobian_plus, vector< vector< vector<double> > > &Jacobian_minus, int nt);




#endif /* OP2A_CFD_BC_HPP_ */
