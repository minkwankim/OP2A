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
#include "CFD/include/VariableSamplesWIG.hpp"


void ApplicationOP2A::preprocessing_grid()
{
	grid.readMeshData(problem_setup.mesh_file_name, static_cast<GRID::GridDataType>(problem_setup.mesh_file_type));
	grid.processingGridData(problem_setup.grid_factor, problem_setup.is_axisymmetric, use_extended_Stencil);
	check_elapsed_time("Reading/Processing Grid data");


	/*
	 * create data template for cells/faces
	 */
	int nBasic_cell_1D	= 10;
	int nBasic_cell_2D	= 0;

	int nBasic_face_1D	= 1;
	int nBasic_face_2D	= 0;

	if (problem_setup.is_axisymmetric != true)	nBasic_cell_1D++;


	if (problem_setup.TIME_INTEGRATION_METHOD != 0)
	{
		nBasic_cell_1D++;
		nBasic_cell_2D++;

		nBasic_face_1D++;
		nBasic_face_2D	+= 3;
	}

	if (problem_setup.is_viscous == true)
	{
		nBasic_cell_1D += 2;
		nBasic_cell_2D++;

		nBasic_face_1D++;
		if (problem_setup.TIME_INTEGRATION_METHOD != 0) nBasic_face_2D	+= 2;
	}

	cell_data1D_template.resize(nBasic_cell_1D);
	cell_data2D_template.resize(nBasic_cell_2D);
	face_data1D_template.resize(nBasic_face_1D);
	face_data2D_template.resize(nBasic_face_2D);



	cell_data1D_template.data[0]	= data_CFD_Q;
	cell_data1D_template.data[1]	= data_CFD_V;
	cell_data1D_template.data[2]	= data_CFD_W;
	cell_data1D_template.data[3]	= data_CFD_mixture;
	cell_data1D_template.data[4]	= data_CFD_Xs;
	cell_data1D_template.data[5]	= data_CFD_Ys;
	cell_data1D_template.data[6]	= data_CFD_R;
	cell_data1D_template.data[7]	= data_CFD_dQ;
	cell_data1D_template.data[8]	= data_CFD_Q;	cell_data1D_template.data[8].asgName("Q_new");
	cell_data1D_template.data[9]	= data_CFD_Source;

	face_data1D_template.data[0]	= data_CFD_Flux_inviscid;


	nBasic_cell_1D	= 10;
	nBasic_cell_2D	= 0;
	nBasic_face_1D	= 1;
	nBasic_face_2D	= 0;

	if (problem_setup.is_axisymmetric != true)
	{
		cell_data1D_template.data[nBasic_cell_1D]	= data_CFD_divVc;	nBasic_cell_1D++;
	}


	if (problem_setup.TIME_INTEGRATION_METHOD != 0)
	{
		cell_data1D_template.data[nBasic_cell_1D]	= data_CFD_dp_dQ;	nBasic_cell_1D++;
		cell_data2D_template.data[nBasic_cell_2D]	= data_CFD_dT_dQ;	nBasic_cell_2D++;

		face_data2D_template.data[nBasic_face_2D]	= data_CFD_Jacobian_inviscid_minus; nBasic_face_2D++;
		face_data2D_template.data[nBasic_face_2D]	= data_CFD_Jacobian_inviscid_plus; nBasic_face_2D++;
		face_data2D_template.data[nBasic_face_2D]	= data_CFD_Jacobian_source; nBasic_face_2D++;
	}


	if (problem_setup.is_viscous == true)
	{
		cell_data1D_template.data[nBasic_cell_1D]	= data_CFD_diffusion_coeff;	nBasic_cell_1D++;
		cell_data1D_template.data[nBasic_cell_1D]	= data_CFD_viscosity_coeff;	nBasic_cell_1D++;

		cell_data2D_template.data[nBasic_cell_2D]	= data_CFD_thermal_conductivity_coeff; nBasic_cell_2D++;

		face_data1D_template.data[nBasic_face_1D]	= data_CFD_Flux_viscous; nBasic_face_1D++;

		if (problem_setup.TIME_INTEGRATION_METHOD != 0)
		{
			face_data2D_template.data[nBasic_face_2D]	= data_CFD_Jacobian_viscous_minus; nBasic_face_2D++;
			face_data2D_template.data[nBasic_face_2D]	= data_CFD_Jacobian_viscous_plus; nBasic_face_2D++;
		}
	}

	cell_data1D_template.mapping();
	cell_data2D_template.mapping();

	face_data1D_template.mapping();
	face_data2D_template.mapping();


	/*
	 * Assign data type for
	 * 	- cell
	 * 	- face
	 * 	- ghost cell
	 */
	// Cell
	Data::DataStorageVector<Data::DataStorage>		cell_data1D_template1;
	Data::DataStorageVector<Data::DataStorage2D>	cell_data2D_template1;

	cell_data1D_template1 = CFD::CFD_DataTemplateWIG::CellData1D(species_set, grid.ND, problem_setup.NER,  problem_setup.NEV,  problem_setup.NEE,  problem_setup.is_viscous,  problem_setup.is_axisymmetric, problem_setup.TIME_INTEGRATION_METHOD);
	cell_data2D_template1 = CFD::CFD_DataTemplateWIG::CellData2D(species_set, grid.ND, problem_setup.NER,  problem_setup.NEV,  problem_setup.NEE,  problem_setup.is_viscous,  problem_setup.is_axisymmetric, problem_setup.TIME_INTEGRATION_METHOD);

	for (int c = 0; c <= grid.NCM; c++)	grid.cells[c].data1D = cell_data1D_template1;
	for (int c = 0; c <= grid.NCM; c++)	grid.cells[c].data2D = cell_data2D_template1;

	// Face
	for (int f = 0; f <= grid.NFM; f++)	grid.faces[f].data1D	= face_data1D_template;
	for (int f = 0; f <= grid.NFM; f++)	grid.faces[f].data2D	= face_data2D_template;




	//line_finder(&grid, grid.grid_line.lines, grid.grid_line.lines_bd, grid.grid_line.cell_line_info, grid.grid_line.num_lines);
}
