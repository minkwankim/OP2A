/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 12, 2015
 *      			Author: Minkwan Kim
 *
 * preprocessing_grid.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_grid.hpp"


/*
 * ================================================================================
 * 		Mesh Read / Processing
 * ================================================================================
 */

void preprocessing_grid_read_process_ver1(string mesh_file_name, int grid_type, bool is_axisymmetric, double grid_factor, GRID_CLASS &grid)
{
	vector<int>	bc_zone(MAX_ZONE_READ_MESH);

	read_mesh_info(mesh_file_name, grid.DIM, grid.NNM, grid.NFM ,grid.NCM, grid.NGM, bc_zone, grid_type);
	read_mesh_node(mesh_file_name, grid.DIM, grid.NNM, grid.nodes.data_ptr, grid.nodes.whereis, grid_type);
	read_mesh_face(mesh_file_name, grid.DIM, grid.NFM, grid.faces.data_ptr, grid.faces.whereis, bc_zone, grid_type);
	read_mesh_cell(mesh_file_name, grid.DIM, grid.NCM, grid.cells.data_ptr, grid.cells.whereis, grid_type);

	processing_node_data(grid.DIM, grid.NNM, grid.nodes.data_ptr, grid_factor, is_axisymmetric);
	processing_face_data(grid.DIM, grid.NFM, grid.nodes.data_ptr, grid.nodes.whereis, grid.faces.data_ptr);
	processing_cell_data(grid.DIM, grid.NFM, grid.NCM, grid.nodes.data_ptr, grid.nodes.whereis, grid.faces.data_ptr, grid.faces.whereis, grid.cells.data_ptr, grid.cells.whereis);
	processing_ghost_cell_data(grid.DIM, grid.NFM, grid.NCM, grid.NGM, grid.nodes.data_ptr, grid.nodes.whereis, grid.faces.data_ptr, grid.faces.whereis, grid.cells.data_ptr, grid.cells.whereis, grid.cells_ghost.data_ptr, grid.cells_ghost.whereis);


	processing_stencil_ver1(grid.DIM, grid.NNM, grid.NFM, grid.NCM, grid.NGM,
							grid.nodes.data_ptr, grid.nodes.whereis,
							grid.faces.data_ptr, grid.faces.whereis,
							grid.cells.data_ptr, grid.cells.whereis,
							grid.cells_ghost.data_ptr, grid.cells_ghost.whereis);


	processing_cell_data_distance_to_wall(grid.DIM, grid.NFM, grid.NCM, grid.NGM, grid.nodes.data_ptr, grid.nodes.whereis, grid.faces.data_ptr, grid.faces.whereis, grid.cells.data_ptr, grid.cells.whereis, grid.cells_ghost.data_ptr, grid.cells_ghost.whereis);;
}
