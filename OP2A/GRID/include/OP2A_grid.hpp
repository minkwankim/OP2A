/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 9, 2015
 *      			Author: Minkwan Kim
 *
 * OPPA_grid.hpp
 * 			-  
 *  
 */

#ifndef OPPA_GRID_HPP_
#define OPPA_GRID_HPP_

#include "constants_grid.hpp"
#include "grid_class.hpp"
#include "read_grid.hpp"
#include "processing_grid.hpp"


void preprocessing_grid_read_process_ver1(string mesh_file_name, int grid_type, bool is_axisymmetric, double grid_factor, GRID_CLASS &grid);


#define read_mesh_info		read_mesh_info_ver2
#define read_mesh_node		read_mesh_node_ver2
#define read_mesh_face		read_mesh_face_ver2
#define read_mesh_cell		read_mesh_cell_ver2

#define preprocessing_grid__read_process	preprocessing_grid__read_process_ver1

#endif /* OPPA_GRID_HPP_ */
