/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 9, 2015
 *      			Author: Minkwan Kim
 *
 * PrintResult.cpp
 * 			-  
 *  
 */

#include <fstream>
#include "GRID/include/PrintResult.hpp"

namespace OP2A{
namespace GRID{

void ResultDataPrintTecplotCell(int P, Grid& grid, const string& title, const string& file_name, string& data_to_print)
{

	ofstream tecplot_file(file_name.c_str());


	// FILE HEADER //
	string		filetype;
	filetype	= "FULL";

	tecplot_file << "TITLE = \""	<< title		<< "\"" << endl;
	tecplot_file << "FILETYPE = "	<< filetype				<< endl;

	tecplot_file << "VARIABLES = ";
	tecplot_file << "\"x-direction, [m]\"  \"y-direction, [m]\" ";
	if (grid.ND == 3)	tecplot_file << "\"z-direction, [m]\" ";

	int data_index = grid.cells[1].data1D.dataMap.find(data_to_print);
	int NV	= grid.cells[1].data1D.data[data_index].numData;

	for (int v = 0; v <= NV-1; v++)
	{
		tecplot_file << "\""<< grid.cells[1].data1D.data[data_index].dataMap.getKey(v) << " \"  ";
	}
	tecplot_file << endl;


	// ZONE RECORD //
	tecplot_file << "ZONE"					<< endl;
	tecplot_file << "T = \""				<< "CPU number: " << P	<< "\"" << endl;	// ZONE TITLE
	tecplot_file << "NODES = "				<< grid.NNM	<< endl;						// NUMBER OF NODES
	tecplot_file << "ELEMENTS = "			<< grid.NCM	<< endl;						// NUMBER OF ELEMENTS
	tecplot_file << "DATAPACKING = BLOCK"	<< endl;

	if (grid.ND == 2)	tecplot_file << "ZONETYPE = FEQUADRILATERAL"	<< endl;
	else				tecplot_file << "ZONETYPE = FEBRICK"			<< endl;

	int NVAR	= 0;

	if (NV == 1)
	{
		tecplot_file << "VARLOCATION =([" << 1 + grid.ND	<<"]=CELLCENTERED)" << endl;
	}
	else
	{
		tecplot_file << "VARLOCATION =([" << grid.ND+1 << "-" << NV + grid.ND <<"]=CELLCENTERED)" << endl;
	}


	// WRITE DATA //
	// 1. NODE DATA
	// Write temporary files
	for (int k = 0; k <= grid.ND-1; k++)
	{
		for (int n = 1; n <= grid.NNM; n++)
		{
			tecplot_file << grid.nodes[n].geo.x[k] << endl;
		}

		tecplot_file << endl;
	}


	// 2. FLOW DATA
	for (int v = 0; v <= NV-1; v++)
	{
		for (int c = 1; c <= grid.NCM; c++)
		{
			tecplot_file << grid.cells[c].data1D.data[data_index].data[v] <<endl;
		}
	}



	// 3. WRITE CONNECTIVITY
	if (grid.ND == 2)
	{
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
	}
	else
	{
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
	} // END writing connectivity

	tecplot_file.close();
}



}
}

