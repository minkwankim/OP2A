/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 9, 2015
 *      			Author: Minkwan Kim
 *
 * processing_grid.hpp
 * 			-  
 *  
 */

#ifndef PROCESSING_GRID_HPP_
#define PROCESSING_GRID_HPP_

#include "../include/constants_grid.hpp"
#include "../include/node_class.hpp"
#include "../include/face_class.hpp"
#include "../include/cell_class.hpp"
#include "../include/grid_class.hpp"

void processing_node_data(unsigned int DIM, unsigned int NNM, vector<NODE_CLASS *>	&nodes, double mesh_factor, bool is_axisymmetric);
void processing_face_data(unsigned int DIM, unsigned int NFM, vector<NODE_CLASS *>	&nodes, vector<int> &whereis_node, vector<FACE_CLASS *>	&faces);
void processing_cell_data(unsigned int DIM, unsigned int NFM, unsigned int NCM, vector<NODE_CLASS *>	&nodes, vector<int> &whereis_node, vector<FACE_CLASS *>	&faces, vector<int> &whereis_face, vector<CELL_CLASS *>	&cells, vector<int> &whereis_cell);
void processing_ghost_cell_data_ver1(unsigned int DIM, unsigned int NFM, unsigned int NCM, unsigned int NGM, vector<NODE_CLASS *>	&nodes, vector<int> &whereis_node, vector<FACE_CLASS *>	&faces, vector<int> &whereis_face, vector<CELL_CLASS *>	&cells, vector<int> &whereis_cell, vector<CELL_CLASS *>	&cells_ghost, vector<int> &whereis_cell_ghost);
void processing_stencil_ver1(unsigned int DIM, unsigned int NNM, unsigned int NFM, unsigned int NCM, unsigned int NGM, vector<NODE_CLASS *>	&nodes, vector<int> &whereis_node, vector<FACE_CLASS *>	&faces, vector<int> &whereis_face, vector<CELL_CLASS *>	&cells, vector<int> &whereis_cell, vector<CELL_CLASS *>	&cells_ghost, vector<int> &whereis_cell_ghost);
void processing_cell_data_distance_to_wall(unsigned int DIM, unsigned int NFM, unsigned int NCM, unsigned int NGM, vector<NODE_CLASS *>	&nodes, vector<int> &whereis_node, vector<FACE_CLASS *>	&faces, vector<int> &whereis_face, vector<CELL_CLASS *>	&cells, vector<int> &whereis_cell, vector<CELL_CLASS *>	&cells_ghost, vector<int> &whereis_cell_ghost);

void processing_face_data2(unsigned int DIM, unsigned int NFM, vector<NODE_CLASS>	&nodes, vector<FACE_CLASS>	&faces);
void processing_cell_data2(unsigned int DIM, unsigned int NFM, unsigned int NCM, vector<NODE_CLASS>	&nodes, vector<FACE_CLASS>	&faces, vector<CELL_CLASS>	&cells);



bool determine_line_starting_point(unsigned int BC);
void line_finder(GRID_CLASS *grid, vector <vector <int> > &lines, vector <vector <int> > &lines_bd, vector <int> &cell_line_ID, int &line_num);


#define processing_ghost_cell_data	processing_ghost_cell_data_ver1
#define processing_stencil			processing_stencil_ver1

#endif /* PROCESSING_GRID_HPP_ */
