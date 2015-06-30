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




#include "CFD/include/VariableConstants.hpp"

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

	// Ghost cells
	for (int c = 0; c <= grid.NGM; c++)	grid.cells_ghost[c].data1D = cell_data1D_template;
	for (int c = 0; c <= grid.NGM; c++)	grid.cells_ghost[c].data2D = cell_data2D_template;


	//line_finder(&grid, grid.grid_line.lines, grid.grid_line.lines_bd, grid.grid_line.cell_line_info, grid.grid_line.num_lines);
	indexQ 		= grid.cells[1].data1D.dataMap.find(NAME_Q);
	indexV 		= grid.cells[1].data1D.dataMap.find(NAME_V);
	indexW 		= grid.cells[1].data1D.dataMap.find(NAME_W);
	indexMIX 	= grid.cells[1].data1D.dataMap.find(NAME_MIX);
	indexXs 	= grid.cells[1].data1D.dataMap.find(NAME_XS);
	indexYs 	= grid.cells[1].data1D.dataMap.find(NAME_YS);
	indexResidue  = grid.cells[1].data1D.dataMap.find(NAME_R);
	indexdQ  	= grid.cells[1].data1D.dataMap.find(NAME_dQ);
	//indexQnew  	= grid.cells[1].data1D.dataMap.find(NAME_Qnew);




	check_elapsed_time("Reading/Processing Grid data");
}
