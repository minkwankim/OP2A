/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 23, 2015
 *      			Author: Minkwan Kim
 *
 * ResultDataPrintTecplotNode.cpp
 * 			-  
 *  
 */



#include <fstream>
#include "GRID/include/PrintResult.hpp"

namespace OP2A{
namespace GRID{

void ResultDataPrintTecplotNode2D(int P, Grid& grid, const string& title, const string& file_name, const string& data_to_print)
{
	ofstream tecplot_file(file_name.c_str());
	int data_index = grid.cells[1].data1D.dataMap.find(data_to_print);
	int NV	= grid.cells[1].data1D.data[data_index].numData;


	// Calculating node properties
	vector<Data::DataStorage>	NodeData(grid.NNM+1, grid.cells[1].data1D.data[data_index]);

#pragma omp parallel for
	for (int n = 1; n <= grid.NNM; n++)
	{
		for (int v = 0; v <= NV-1; v++)
		{
			NodeData[n].data[v] = 0.0;
			int numCell = grid.nodes[n].geo.NSC;
			for (int c = 0; c <= numCell-1; c++)
			{
				NodeData[n].data[v] += grid.nodes[n].geo.CellList[c]->data1D.data[data_index].data[v] * grid.nodes[n].geo.WeightShairedCell[c];
			}
		}
	}



	// FILE HEADER //
	string		filetype;
	filetype	= "FULL";

	tecplot_file << "TITLE = \""	<< title		<< "\"" << endl;
	tecplot_file << "FILETYPE = "	<< filetype				<< endl;

	tecplot_file << "VARIABLES = ";
	tecplot_file << "\"x-direction, [m]\"  \"y-direction, [m]\" ";
	for (int v = 0; v <= NV-1; v++)
	{
		tecplot_file << "\""<< grid.cells[1].data1D.data[data_index].dataMap.findKey(v) << " \"  ";
	}
	tecplot_file << endl;


	// ZONE RECORD //
	tecplot_file << "ZONE"					<< endl;
	tecplot_file << "T = \""				<< "CPU number: " << P	<< "\"" << endl;	// ZONE TITLE
	tecplot_file << "NODES = "				<< grid.NNM	<< endl;						// NUMBER OF NODES
	tecplot_file << "ELEMENTS = "			<< grid.NCM	<< endl;						// NUMBER OF ELEMENTS
	tecplot_file << "F = FEPOINT"	<< endl;
	tecplot_file << "ET = QUADRILATERAL"	<< endl;


	// WRITE DATA //
	for (int n = 1; n <= grid.NNM; n++)
	{
		for (int k = 0; k <= grid.ND-1; k++)	tecplot_file << grid.nodes[n].geo.x[k] << " ";

		for (int v = 0; v <= NV-1; v++)			tecplot_file << NodeData[n].data[v] << " ";
		tecplot_file << endl;
	}

	// WRITE CONNECTIVITY
	for (int c = 1; c <= grid.NCM; c++)
	{
		switch (grid.cells[c].geo.type)
		{
		case CellType::quadrilateral:
			for (int n = 0; n <= 3; n++)
			{
				tecplot_file << grid.cells[c].geo.node_list[n]->geo.ID << " ";
			}
			break;

		case CellType::triangle:
			for (int n = 0; n <= 2; n++)
			{
				tecplot_file << grid.cells[c].geo.node_list[n]->geo.ID << " ";
			}
			tecplot_file << grid.cells[c].geo.node_list[2]->geo.ID << " ";
			break;
		}

		tecplot_file << endl;
	}

	tecplot_file.close();
}





void ResultDataPrintTecplotNode3D(int P, Grid& grid, const string& title, const string& file_name, const string& data_to_print)
{
	ofstream tecplot_file(file_name.c_str());
	int data_index = grid.cells[1].data1D.dataMap.find(data_to_print);
	int NV	= grid.cells[1].data1D.data[data_index].numData;


	// Calculating node properties
	vector<Data::DataStorage>	NodeData(grid.NNM+1, grid.cells[1].data1D.data[data_index]);

#pragma omp parallel for
	for (int n = 1; n <= grid.NNM; n++)
	{
		for (int v = 0; v <= NV-1; v++)
		{
			NodeData[n].data[v] = 0.0;
			int numCell = grid.nodes[n].geo.NSC;
			for (int c = 0; c <= numCell-1; c++)
			{
				NodeData[n].data[v] += grid.nodes[n].geo.CellList[c]->data1D.data[data_index].data[v] * grid.nodes[n].geo.WeightShairedCell[c];
			}
		}
	}



	// FILE HEADER //
	string		filetype;
	filetype	= "FULL";

	tecplot_file << "TITLE = \""	<< title		<< "\"" << endl;
	tecplot_file << "FILETYPE = "	<< filetype				<< endl;

	tecplot_file << "VARIABLES = ";
	tecplot_file << "\"x-direction, [m]\"  \"y-direction, [m]\"  \"z-direction, [m]\" ";
	for (int v = 0; v <= NV-1; v++)
	{
		tecplot_file << "\""<< grid.cells[1].data1D.data[data_index].dataMap.findKey(v) << " \"  ";
	}
	tecplot_file << endl;


	// Internal Flow
	// 	- ZONE RECORD: Internal Flow//
	tecplot_file << "ZONE"		 	<< endl;
	tecplot_file << "T = FLUID "	 << endl;	// ZONE TITLE
	tecplot_file << "N = "			<< grid.NNM	<< endl;						// NUMBER OF NODES
	tecplot_file << "E = "			<< grid.NCM	<< endl;						// NUMBER OF ELEMENTS
	tecplot_file << "F = FEPOINT"	<< endl;
	tecplot_file << "ET = BRICK"	<< endl;

	// 	- WRITE DATA //
	for (int n = 1; n <= grid.NNM; n++)
	{
		for (int k = 0; k <= grid.ND-1; k++)	tecplot_file << grid.nodes[n].geo.x[k] << " ";
		for (int v = 0; v <= NV-1; v++)			tecplot_file << NodeData[n].data[v] << " ";
		tecplot_file << endl;
	}

	//  - WRITE CONNECTIVITY
	for (int c = 1; c <= grid.NCM; c++)
	{
		switch (grid.cells[c].geo.type)
		{
		case CellType::tetrahedron:
			tecplot_file << grid.cells[c].geo.node_list[0]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[1]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[2]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[2]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[3]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[3]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[3]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[3]->geo.ID;
			tecplot_file << endl;
			break;

		case CellType::hexahedron:
			tecplot_file << grid.cells[c].geo.node_list[0]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[1]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[2]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[3]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[4]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[5]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[6]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[7]->geo.ID;
			tecplot_file << endl;
			break;

		case CellType::pyramid:
			tecplot_file << grid.cells[c].geo.node_list[0]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[1]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[2]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[3]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[4]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[4]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[4]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[4]->geo.ID;
			tecplot_file << endl;
			break;

		case CellType::wedge:
			tecplot_file << grid.cells[c].geo.node_list[0]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[1]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[2]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[2]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[3]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[4]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[5]->geo.ID <<" ";
			tecplot_file << grid.cells[c].geo.node_list[6]->geo.ID;
			tecplot_file << endl;
			break;
		}
	}


	int numFacesWall 			= 0;
	int numFacesInlet 			= 0;
	int numFacesOutlet 			= 0;
	int numFacesFreestream 		= 0;
	int numFacesSymmetric 		= 0;
	int numFacesAxis 			= 0;
	int numFacesAnode 			= 0;
	int numFacesCathode 		= 0;
	int numFacesDielectricWall 	= 0;

	for (int c = 1; c <= grid.NGM; c++)
	{
		switch (grid.cells_ghost[c].geo.BC)
		{
		case BCType::wall:
			numFacesWall++;
			break;
		case BCType::inlet:
			numFacesInlet++;
			break;
		case BCType::outlet:
			numFacesOutlet++;
			break;
		case BCType::freestream:
			numFacesFreestream++;
			break;
		case BCType::symmetric:
			numFacesSymmetric++;
			break;
		case BCType::axis:
			numFacesAxis++;
			break;
		case BCType::anode:
			numFacesAnode++;
			break;
		case BCType::cathode:
			numFacesCathode++;
			break;
		case BCType::dielectricwall:
			numFacesDielectricWall++;
			break;
		}
	}

	// 1. Wall
	if (numFacesWall > 0)
	{
		tecplot_file << "ZONE"		 			<< endl;
		tecplot_file << "T = BODY "	 			<< endl;	// ZONE TITLE
		tecplot_file << "N = "					<< grid.NNM		<< endl;						// NUMBER OF NODES
		tecplot_file << "E = "					<< numFacesWall	<< endl;						// NUMBER OF ELEMENTS
		tecplot_file << "F = FEPOINT"			<< endl;
		tecplot_file << "ET = QUADRILATERAL"	<< endl;

		// 	- WRITE DATA //
		for (int n = 1; n <= grid.NNM; n++)
		{
			for (int k = 0; k <= grid.ND-1; k++)	tecplot_file << grid.nodes[n].geo.x[k] << " ";
			for (int v = 0; v <= NV-1; v++)			tecplot_file << NodeData[n].data[v] << " ";
			tecplot_file << endl;

			for (int c = 1; c <= grid.NGM; c++)
			{
				if (grid.cells_ghost[c].geo.BC == int(BCType::wall))
				{
					switch (grid.cells_ghost[c].geo.face_list[0]->geo.type)
					{
					case FaceType::f_triangle:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					case FaceType::f_quadrilateral:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[3]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					}
				}
			}
		}
	}


	// 2. Inlet
	if (numFacesInlet > 0)
	{
		tecplot_file << "ZONE"		 			<< endl;
		tecplot_file << "T = INLET "	 		<< endl;	// ZONE TITLE
		tecplot_file << "N = "					<< grid.NNM		<< endl;						// NUMBER OF NODES
		tecplot_file << "E = "					<< numFacesInlet	<< endl;						// NUMBER OF ELEMENTS
		tecplot_file << "F = FEPOINT"			<< endl;
		tecplot_file << "ET = QUADRILATERAL"	<< endl;

		// 	- WRITE DATA //
		for (int n = 1; n <= grid.NNM; n++)
		{
			for (int k = 0; k <= grid.ND-1; k++)	tecplot_file << grid.nodes[n].geo.x[k] << " ";
			for (int v = 0; v <= NV-1; v++)			tecplot_file << NodeData[n].data[v] << " ";
			tecplot_file << endl;

			for (int c = 1; c <= grid.NGM; c++)
			{
				if (grid.cells_ghost[c].geo.BC == int(BCType::inlet))
				{
					switch (grid.cells_ghost[c].geo.face_list[0]->geo.type)
					{
					case FaceType::f_triangle:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					case FaceType::f_quadrilateral:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[3]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					}
				}
			}
		}
	}


	// 3. Outlet
	if (numFacesOutlet > 0)
	{
		tecplot_file << "ZONE"		 			<< endl;
		tecplot_file << "T = OUTLET "	 		<< endl;	// ZONE TITLE
		tecplot_file << "N = "					<< grid.NNM		<< endl;						// NUMBER OF NODES
		tecplot_file << "E = "					<< numFacesOutlet	<< endl;						// NUMBER OF ELEMENTS
		tecplot_file << "F = FEPOINT"			<< endl;
		tecplot_file << "ET = QUADRILATERAL"	<< endl;

		// 	- WRITE DATA //
		for (int n = 1; n <= grid.NNM; n++)
		{
			for (int k = 0; k <= grid.ND-1; k++)	tecplot_file << grid.nodes[n].geo.x[k] << " ";
			for (int v = 0; v <= NV-1; v++)			tecplot_file << NodeData[n].data[v] << " ";
			tecplot_file << endl;

			for (int c = 1; c <= grid.NGM; c++)
			{
				if (grid.cells_ghost[c].geo.BC == int(BCType::outlet))
				{
					switch (grid.cells_ghost[c].geo.face_list[0]->geo.type)
					{
					case FaceType::f_triangle:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					case FaceType::f_quadrilateral:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[3]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					}
				}
			}
		}
	}


	// 4. Freestream
	if (numFacesFreestream > 0)
	{
		tecplot_file << "ZONE"		 			<< endl;
		tecplot_file << "T = FREESTREAM "	 	<< endl;	// ZONE TITLE
		tecplot_file << "N = "					<< grid.NNM				<< endl;						// NUMBER OF NODES
		tecplot_file << "E = "					<< numFacesFreestream	<< endl;						// NUMBER OF ELEMENTS
		tecplot_file << "F = FEPOINT"			<< endl;
		tecplot_file << "ET = QUADRILATERAL"	<< endl;

		// 	- WRITE DATA //
		for (int n = 1; n <= grid.NNM; n++)
		{
			for (int k = 0; k <= grid.ND-1; k++)	tecplot_file << grid.nodes[n].geo.x[k] << " ";
			for (int v = 0; v <= NV-1; v++)			tecplot_file << NodeData[n].data[v] << " ";
			tecplot_file << endl;

			for (int c = 1; c <= grid.NGM; c++)
			{
				if (grid.cells_ghost[c].geo.BC == int(BCType::freestream))
				{
					switch (grid.cells_ghost[c].geo.face_list[0]->geo.type)
					{
					case FaceType::f_triangle:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					case FaceType::f_quadrilateral:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[3]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					}
				}
			}
		}
	}



	// 5. Symmetric
	if (numFacesSymmetric > 0)
	{
		tecplot_file << "ZONE"		 			<< endl;
		tecplot_file << "T = SYMMETRIC "	 	<< endl;	// ZONE TITLE
		tecplot_file << "N = "					<< grid.NNM				<< endl;						// NUMBER OF NODES
		tecplot_file << "E = "					<< numFacesSymmetric	<< endl;						// NUMBER OF ELEMENTS
		tecplot_file << "F = FEPOINT"			<< endl;
		tecplot_file << "ET = QUADRILATERAL"	<< endl;

		// 	- WRITE DATA //
		for (int n = 1; n <= grid.NNM; n++)
		{
			for (int k = 0; k <= grid.ND-1; k++)	tecplot_file << grid.nodes[n].geo.x[k] << " ";
			for (int v = 0; v <= NV-1; v++)			tecplot_file << NodeData[n].data[v] << " ";
			tecplot_file << endl;

			for (int c = 1; c <= grid.NGM; c++)
			{
				if (grid.cells_ghost[c].geo.BC == int(BCType::symmetric))
				{
					switch (grid.cells_ghost[c].geo.face_list[0]->geo.type)
					{
					case FaceType::f_triangle:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					case FaceType::f_quadrilateral:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[3]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					}
				}
			}
		}
	}


	// 6. Axis
	if (numFacesAxis > 0)
	{
		tecplot_file << "ZONE"		 			<< endl;
		tecplot_file << "T = AXIS "	 		<< endl;	// ZONE TITLE
		tecplot_file << "N = "					<< grid.NNM		<< endl;						// NUMBER OF NODES
		tecplot_file << "E = "					<< numFacesAxis	<< endl;						// NUMBER OF ELEMENTS
		tecplot_file << "F = FEPOINT"			<< endl;
		tecplot_file << "ET = QUADRILATERAL"	<< endl;

		// 	- WRITE DATA //
		for (int n = 1; n <= grid.NNM; n++)
		{
			for (int k = 0; k <= grid.ND-1; k++)	tecplot_file << grid.nodes[n].geo.x[k] << " ";
			for (int v = 0; v <= NV-1; v++)			tecplot_file << NodeData[n].data[v] << " ";
			tecplot_file << endl;

			for (int c = 1; c <= grid.NGM; c++)
			{
				if (grid.cells_ghost[c].geo.BC == int(BCType::axis))
				{
					switch (grid.cells_ghost[c].geo.face_list[0]->geo.type)
					{
					case FaceType::f_triangle:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					case FaceType::f_quadrilateral:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[3]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					}
				}
			}
		}
	}


	// 7. Anode
	if (numFacesAnode > 0)
	{
		tecplot_file << "ZONE"		 			<< endl;
		tecplot_file << "T = ANODE "	 		<< endl;	// ZONE TITLE
		tecplot_file << "N = "					<< grid.NNM		<< endl;						// NUMBER OF NODES
		tecplot_file << "E = "					<< numFacesAnode	<< endl;						// NUMBER OF ELEMENTS
		tecplot_file << "F = FEPOINT"			<< endl;
		tecplot_file << "ET = QUADRILATERAL"	<< endl;

		// 	- WRITE DATA //
		for (int n = 1; n <= grid.NNM; n++)
		{
			for (int k = 0; k <= grid.ND-1; k++)	tecplot_file << grid.nodes[n].geo.x[k] << " ";
			for (int v = 0; v <= NV-1; v++)			tecplot_file << NodeData[n].data[v] << " ";
			tecplot_file << endl;

			for (int c = 1; c <= grid.NGM; c++)
			{
				if (grid.cells_ghost[c].geo.BC == int(BCType::anode))
				{
					switch (grid.cells_ghost[c].geo.face_list[0]->geo.type)
					{
					case FaceType::f_triangle:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					case FaceType::f_quadrilateral:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[3]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					}
				}
			}
		}
	}



	// 8. Cathode
	if (numFacesCathode > 0)
	{
		tecplot_file << "ZONE"		 			<< endl;
		tecplot_file << "T = CATHODE "	 		<< endl;	// ZONE TITLE
		tecplot_file << "N = "					<< grid.NNM			<< endl;						// NUMBER OF NODES
		tecplot_file << "E = "					<< numFacesCathode	<< endl;						// NUMBER OF ELEMENTS
		tecplot_file << "F = FEPOINT"			<< endl;
		tecplot_file << "ET = QUADRILATERAL"	<< endl;

		// 	- WRITE DATA //
		for (int n = 1; n <= grid.NNM; n++)
		{
			for (int k = 0; k <= grid.ND-1; k++)	tecplot_file << grid.nodes[n].geo.x[k] << " ";
			for (int v = 0; v <= NV-1; v++)			tecplot_file << NodeData[n].data[v] << " ";
			tecplot_file << endl;

			for (int c = 1; c <= grid.NGM; c++)
			{
				if (grid.cells_ghost[c].geo.BC == int(BCType::cathode))
				{
					switch (grid.cells_ghost[c].geo.face_list[0]->geo.type)
					{
					case FaceType::f_triangle:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					case FaceType::f_quadrilateral:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[3]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					}
				}
			}
		}
	}



	// 9. Dielectricwall
	if (numFacesDielectricWall > 0)
	{
		tecplot_file << "ZONE"		 			<< endl;
		tecplot_file << "T = DIELECTRICWALL "	 		<< endl;	// ZONE TITLE
		tecplot_file << "N = "					<< grid.NNM		<< endl;						// NUMBER OF NODES
		tecplot_file << "E = "					<< numFacesDielectricWall	<< endl;						// NUMBER OF ELEMENTS
		tecplot_file << "F = FEPOINT"			<< endl;
		tecplot_file << "ET = QUADRILATERAL"	<< endl;

		// 	- WRITE DATA //
		for (int n = 1; n <= grid.NNM; n++)
		{
			for (int k = 0; k <= grid.ND-1; k++)	tecplot_file << grid.nodes[n].geo.x[k] << " ";
			for (int v = 0; v <= NV-1; v++)			tecplot_file << NodeData[n].data[v] << " ";
			tecplot_file << endl;

			for (int c = 1; c <= grid.NGM; c++)
			{
				if (grid.cells_ghost[c].geo.BC == int(BCType::dielectricwall))
				{
					switch (grid.cells_ghost[c].geo.face_list[0]->geo.type)
					{
					case FaceType::f_triangle:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					case FaceType::f_quadrilateral:
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[0]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[1]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[2]->geo.ID <<" ";
						tecplot_file << grid.cells_ghost[c].geo.face_list[0]->geo.node_list[3]->geo.ID <<" ";
						tecplot_file << endl;
						break;
					}
				}
			}
		}
	}


	tecplot_file.close();
}




void ResultDataPrintTecplotNode(int P, Grid& grid, const string& title, const string& file_name, const string& data_to_print)
{
	if (grid.ND == 2)
	{
		ResultDataPrintTecplotNode2D(P, grid, title, file_name, data_to_print);
	}
	else
	{
		ResultDataPrintTecplotNode3D(P, grid, title, file_name, data_to_print);
	}
}




}
}
