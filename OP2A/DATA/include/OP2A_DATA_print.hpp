/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 12, 2015
 *      			Author: Minkwan Kim
 *
 * data_print.hpp
 * 			-  
 *  
 */

#ifndef DATA_PRINT_HPP_
#define DATA_PRINT_HPP_



#include "../include/data_class.hpp"
#include "../../GRID/include/OP2A_grid.hpp"


#define PRINT_NORMAL	0
#define PRINT_TIME		1

#define	TECPLOT	0


/*
 * ==============================================================
 * 		Data print function for Tecplot
 * ==============================================================
 */
void OP2A_data_print_tecplot_ver1(int P, GRID_CLASS &grid_data, vector < vector<double> > &Vc, vector<string> &variable_names, string title, string file_name);
void OP2A_data_print_tecplot_time_ver1(int P, GRID_CLASS &grid_data, vector < vector<double> > &Vc, vector<string> &variable_names, string title, string file_name, double t, int flag);
void OP2A_data_print_tecplot_multi_ver1(int P, GRID_CLASS &grid_data, vector< vector < vector<double> > > &Vc, vector< vector<string> > &variable_names, string title, string file_name, int num_fluid);
void OP2A_data_print_restart(int P, GRID_CLASS &grid_data, vector< vector < vector<double> > > &Vc, string file_name, int num_fluid, int iter, int type);

#define OP2A_data_print_tecplot			OP2A_data_print_tecplot_ver1
#define OP2A_data_print_tecplot_time 	OP2A_data_print_tecplot_time_ver1
#define OP2A_data_print_tecplot_multi	OP2A_data_print_tecplot_multi_ver1



/*
 * ==============================================================
 * 		Data print function
 * ==============================================================
 */



#endif /* DATA_PRINT_HPP_ */
