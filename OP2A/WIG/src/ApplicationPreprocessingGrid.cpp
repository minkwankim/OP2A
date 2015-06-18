/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 12, 2015
 *      			Author: Minkwan Kim
 *
 * ApplicationPreprocessingGrid.cpp
 * 			-  
 *  
 */





#include "Common/include/CodeLocation.hpp"
#include "Common/include/Exception_NPExceed.hpp"
#include "Common/include/StringOps.hpp"

#include "../include/OP2A_Application.hpp"


void ApplicationOP2A::preprocessing_grid()
{
	grid.readMeshData(problem_setup.mesh_file_name, static_cast<GRID::GridDataType>(problem_setup.mesh_file_type));
	grid.processingGridData(problem_setup.grid_factor, problem_setup.is_axisymmetric, use_extended_Stencil);


	/*
	 * assign data template for cells/faces
	 */
	// Cell
	for (int c = 0; c <= grid.NCM; c++)	grid.cells[c].data1D = cell_data1D_template;
	for (int c = 0; c <= grid.NCM; c++)	grid.cells[c].data2D = cell_data2D_template;

	// Face
	for (int f = 0; f <= grid.NFM; f++)	grid.faces[f].data1D	= face_data1D_template;
	for (int f = 0; f <= grid.NFM; f++)	grid.faces[f].data2D	= face_data2D_template;


	//line_finder(&grid, grid.grid_line.lines, grid.grid_line.lines_bd, grid.grid_line.cell_line_info, grid.grid_line.num_lines);
	check_elapsed_time("Reading/Processing Grid data");
}
