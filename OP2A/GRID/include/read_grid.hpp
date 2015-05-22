/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 9, 2015
 *      			Author: Minkwan Kim
 *
 * read_grid.cpp
 * 			-  
 *  
 */

#ifndef READ_GRID_CPP_
#define READ_GRID_CPP_

#include <vector>
#include <string>
using namespace std;





/*
 * ================================================================================
 * 		READ Basic grid information from grid file
 * ================================================================================
 */

void read_mesh_info_ver2(string mesh_file_name,		// Mesh file name
						unsigned int& DIM,			// Dimention
						unsigned int& NNM,			// Total number of nodes in mesh
					 	unsigned int& NFM,			// Total number of faces in mesh
					 	unsigned int& NCM,			// Total number of cells in mesh
					 	unsigned int& NGM,			// Total number of ghost-cells in mesh
					 	vector<int>& bc_zone,		// BC
					 	int type);

void read_mesh_info_fluent(string mesh_file_name,		// Mesh file name
							unsigned int& DIM,			// Dimenstion
							unsigned int& NNM,			// Total number of nodes in mesh
							unsigned int& NFM,			// Total number of faces in mesh
							unsigned int& NCM,			// Total number of cells in mesh
							unsigned int& NGM,			// Total number of ghost-cells in mesh
							vector<int>& bc_zone);		// BC Information





/*
 * ================================================================================
 * 		READ node data
 * ================================================================================
 */
void read_mesh_node_ver2(string 		mesh_file_name,		// Mesh file name
						unsigned int	DIM,
						unsigned int	NNM,				// Grid information
						vector<NODE_CLASS *>	&nodes,
						vector<int>			&whereis,
						int TYPE);

void read_mesh_node_fluent(string 		mesh_file_name,		// Mesh file name
							unsigned int	DIM,
							unsigned int	NNM,				// Grid information
							vector<NODE_CLASS *>	&nodes,
							vector<int>				&whereis);





/*
 * ================================================================================
 * 		READ face data
 * ================================================================================
 */
void read_mesh_face_ver2(string 		mesh_file_name,		// Mesh file name
						unsigned int	DIM,
						unsigned int	NFM,				// Grid information
						vector<FACE_CLASS *>	&faces,
						vector<int>				&whereis,
						vector<int>				&bc_zone,
						int TYPE);

void read_mesh_face_fluent(string 		mesh_file_name,		// Mesh file name
							unsigned int	DIM,
							unsigned int	NFM,				// Grid information
							vector<FACE_CLASS *>	&faces,
							vector<int>				&whereis,
							vector<int>				&bc_zone);





/*
 * ================================================================================
 * 		READ cell data
 * ================================================================================
 */
void read_mesh_cell_ver2(string mesh_file_name,		// Mesh file name
						unsigned int	DIM,
						unsigned int	NCM,				// Grid information
						vector<CELL_CLASS *>	&cells,
						vector<int>				&whereis,
						int TYPE);

void read_mesh_cell_fluent(string mesh_file_name,		// Mesh file name
		unsigned int	DIM,
		unsigned int	NCM,				// Grid information
		vector<CELL_CLASS *>	&cells,
		vector<int>				&whereis);







#endif /* READ_GRID_CPP_ */
