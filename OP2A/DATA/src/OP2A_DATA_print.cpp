/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 12, 2015
 *      			Author: Minkwan Kim
 *
 * data_print.cpp
 * 			-  
 *  
 */


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include "../include/OP2A_data.hpp"
#include "../../GRID/include/OP2A_grid.hpp"


/*
 * ==============================================================
 * 		Data print function for Tecplot
 * ==============================================================
 */
void OP2A_data_print_tecplot_ver1(int P, GRID_CLASS &grid_data, vector < vector<double> > &Vc, vector<string> &variable_names, string title, string file_name)
{
	int				i, j;
	int				NCM_p = 0;
	int 			node_increase;
	double			var_value;
	string			filetype;
	ofstream		tecplot_file;



	// OPEN FILE TO WRITE
	tecplot_file.open(file_name.c_str());


	// FILE HEADER //
	filetype	= "FULL";
	tecplot_file << "TITLE = \""	<< title		<< "\"" << endl;
	tecplot_file << "FILETYPE = "	<< filetype				<< endl;

	tecplot_file << "VARIABLES = ";
	tecplot_file << "\"x-direction, [m]\"  \"y-direction, [m]\" ";
	if (grid_data.DIM == 3)	tecplot_file << "\"z-direction, [m]\" ";

	int NV = variable_names.size();
	for (i = 0; i <= NV-1; i++)	tecplot_file << "\""<< variable_names[i] << " \"  ";
	tecplot_file << endl;


	// ZONE RECORD //
	tecplot_file << "ZONE"					<< endl;
	tecplot_file << "T = \""				<< "CPU number: " << P	<< "\"" << endl;	// ZONE TITLE
	tecplot_file << "NODES = "				<< grid_data.NNM				<< endl;			// NUMBER OF NODES
	tecplot_file << "ELEMENTS = "			<< grid_data.NCM				<< endl;			// NUMBER OF ELEMENTS
	tecplot_file << "DATAPACKING = BLOCK"	<< endl;

	if (grid_data.DIM == 2)	tecplot_file << "ZONETYPE = FEQUADRILATERAL"	<< endl;
	else					tecplot_file << "ZONETYPE = FEBRICK"			<< endl;

	if (NV == 1)	tecplot_file << "VARLOCATION =([" << NV + grid_data.DIM	<<"]=CELLCENTERED)" << endl;
	else			tecplot_file << "VARLOCATION =([" << grid_data.DIM+1	<< "-" << NV + grid_data.DIM <<"]=CELLCENTERED)" << endl;



	// WRITE DATA //
	// 1. NODE DATA
	// Write temporary files
	int node_ID_ptr;
	for (j = 0; j <= grid_data.DIM-1; j++)
	{
		for (i = 1; i <= grid_data.NNM; i++)
		{
			node_ID_ptr	= grid_data.nodes.whereis[i];
			tecplot_file << grid_data.nodes.data_ptr[node_ID_ptr]->x[j] << endl;
		}

		tecplot_file << endl;
	}


	// 2. FLOW DATA
	int cell_ID_ptr;
	for (j = 0; j <= NV-1; j++)
	{
		for (i = 1; i <= grid_data.NCM; i++)
		{
			cell_ID_ptr	= grid_data.cells.whereis[i];
			tecplot_file << Vc[cell_ID_ptr][j] << endl;
		}
	}


	// 3. WRITE CONNECTIVITY
	if (grid_data.DIM == 2)
	{
		for (i = 1; i <= grid_data.NCM; i++)
		{
			cell_ID_ptr	= grid_data.cells.whereis[i];

			switch (grid_data.cells.data_ptr[cell_ID_ptr]->type)
			{
			case C_QUADRILATERAL:
				for (j = 0; j <= 3; j++)
				{
					tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[j] << " ";
				}
				break;
			case C_TRIANGLE:
				for (j = 0; j <= 2; j++)
				{
					tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[j] << " ";
				}
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2];
				break;
			}

			tecplot_file << endl;
		}
	}
	else
	{
		for (i = 1; i <= grid_data.NCM; i++)
		{
			cell_ID_ptr	= grid_data.cells.whereis[i];
			switch (grid_data.cells.data_ptr[cell_ID_ptr]->type)
			{

			case C_TETRAHEDRON:
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[0] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[1] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3];
				tecplot_file << endl;
				break;

			case C_HEXAHEDRON:
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[0] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[1] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[5] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[6] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[7];
				tecplot_file << endl;
				break;

			case C_PYRAMID:
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[0] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[1] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4];
				tecplot_file << endl;
				break;

			case C_PRISM:
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[0] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[1] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[5] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[6];
				tecplot_file << endl;
				break;
			}
		}
	} // END writing connectivity

	tecplot_file.close();
}

void OP2A_data_print_tecplot_time_ver1(int P, GRID_CLASS &grid_data, vector < vector<double> > &Vc, vector<string> &variable_names, string title, string file_name, double t, int flag)
{
	int		i, j;
	int		NCM_p;
	int 	node_increase;

	double	var_value;

	string			filetype;
	ofstream		tecplot_file;

	NCM_p	= 0;
	int NV = variable_names.size();


	if (flag == 0)
	{
		// OPEN FILE TO WRITE
		tecplot_file.open(file_name.c_str());


		// FILE HEADER //
		filetype	= "FULL";
		tecplot_file << "TITLE = \""	<< title		<< "\"" << endl;
		tecplot_file << "FILETYPE = "	<< filetype				<< endl;

		tecplot_file << "VARIABLES = ";
		tecplot_file << "\"x-direction, [m]\"  \"y-direction, [m]\" ";
		if (grid_data.DIM == 3)	tecplot_file << "\"z-direction, [m]\" ";

		for (i = 0; i <= NV-1; i++)	tecplot_file << "\""<< variable_names[i] << " \"  ";
		tecplot_file << endl;
	}
	else
	{
		tecplot_file.open(file_name.c_str(), ios::app);
		tecplot_file << endl;
	}


	// ZONE RECORD //
	tecplot_file << "ZONE"					<< endl;
	tecplot_file << "T = \""				<< "CPU number: " << P	<< "\"" << endl;	// ZONE TITLE
	tecplot_file << "STRANDID = "			<< flag+1	<< endl;						// TIME STRAND ID
	tecplot_file << "SOLUTIONTIME = "		<< t		<< endl;						// Solution time
	tecplot_file << "NODES = "				<< grid_data.NNM				<< endl;			// NUMBER OF NODES
	tecplot_file << "ELEMENTS = "			<< grid_data.NCM				<< endl;			// NUMBER OF ELEMENTS
	tecplot_file << "DATAPACKING = BLOCK"	<< endl;

	if (grid_data.DIM == 2)	tecplot_file << "ZONETYPE = FEQUADRILATERAL"	<< endl;
	else					tecplot_file << "ZONETYPE = FEBRICK"			<< endl;

	if (NV == 1)	tecplot_file << "VARLOCATION =([" << NV+grid_data.DIM	<<"]=CELLCENTERED)" << endl;
	else			tecplot_file << "VARLOCATION =([" << grid_data.DIM+1	<< "-" << NV+grid_data.DIM <<"]=CELLCENTERED)" << endl;


	// WRITE DATA //
	// 1. NODE DATA
	// Write temporary files
	int node_ID_ptr;
	for (j = 0; j <= grid_data.DIM-1; j++)
	{
		for (i = 1; i <= grid_data.NNM; i++)
		{
			node_ID_ptr	= grid_data.nodes.whereis[i];
			tecplot_file << grid_data.nodes.data_ptr[node_ID_ptr]->x[j] << endl;
		}

		tecplot_file << endl;
	}


	// 2. FLOW DATA
	int cell_ID_ptr;
	for (j = 0; j <= NV-1; j++)
	{
		for (i = 1; i <= grid_data.NCM; i++)
		{
			cell_ID_ptr	= grid_data.cells.whereis[i];
			tecplot_file << Vc[cell_ID_ptr][j] <<endl;
		}
	}


	// 3. WRITE CONNECTIVITY
	if (grid_data.DIM == 2)
	{
		for (i = 1; i <= grid_data.NCM; i++)
		{
			cell_ID_ptr	= grid_data.cells.whereis[i];

			switch (grid_data.cells.data_ptr[cell_ID_ptr]->type)
			{
			case C_QUADRILATERAL:
				for (j = 0; j <= 3; j++)
				{
					tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[j] << " ";
				}
				break;
			case C_TRIANGLE:
				for (j = 0; j <= 2; j++)
				{
					tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[j] << " ";
				}
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2];
				break;
			}

			tecplot_file << endl;
		}
	}
	else
	{
		for (i = 1; i <= grid_data.NCM; i++)
		{
			cell_ID_ptr	= grid_data.cells.whereis[i];
			switch (grid_data.cells.data_ptr[cell_ID_ptr]->type)
			{

			case C_TETRAHEDRON:
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[0] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[1] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3];
				tecplot_file << endl;
				break;

			case C_HEXAHEDRON:
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[0] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[1] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[5] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[6] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[7];
				tecplot_file << endl;
				break;

			case C_PYRAMID:
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[0] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[1] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4];
				tecplot_file << endl;
				break;

			case C_PRISM:
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[0] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[1] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[5] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[6];
				tecplot_file << endl;
				break;
			}
		}
	} // END writing connectivity

	tecplot_file.close();
}



void OP2A_data_print_tecplot_multi_ver1(int P, GRID_CLASS &grid_data, vector< vector < vector<double> > > &Vc, vector< vector<string> > &variable_names, string title, string file_name, int num_fluid)
{
	int		i, j;
	int 	NV;
	int		NCM_p	= 0;
	int 	node_increase;

	double	var_value;

	string			filetype;
	ofstream		tecplot_file;

	// OPEN FILE TO WRITE
	tecplot_file.open(file_name.c_str());


	// FILE HEADER //
	filetype	= "FULL";
	tecplot_file << "TITLE = \""	<< title		<< "\"" << endl;
	tecplot_file << "FILETYPE = "	<< filetype				<< endl;

	tecplot_file << "VARIABLES = ";
	tecplot_file << "\"x-direction, [m]\"  \"y-direction, [m]\" ";
	if (grid_data.DIM == 3)	tecplot_file << "\"z-direction, [m]\" ";

	for (int f = 0; f <= num_fluid-1; f++)
	{
		NV = variable_names[f].size();
		for (i = 0; i <= NV-1; i++)	tecplot_file << "\""<< variable_names[f][i] << "<sub> fluid#" << f << "</sub> \"  ";
	}
	tecplot_file << endl;


	// ZONE RECORD //
	tecplot_file << "ZONE"					<< endl;
	tecplot_file << "T = \""				<< "CPU number: " << P	<< "\"" << endl;	// ZONE TITLE
	tecplot_file << "NODES = "				<< grid_data.NNM				<< endl;			// NUMBER OF NODES
	tecplot_file << "ELEMENTS = "			<< grid_data.NCM				<< endl;			// NUMBER OF ELEMENTS
	tecplot_file << "DATAPACKING = BLOCK"	<< endl;

	if (grid_data.DIM == 2)	tecplot_file << "ZONETYPE = FEQUADRILATERAL"	<< endl;
	else					tecplot_file << "ZONETYPE = FEBRICK"			<< endl;

	int NVAR	= 0;

	if (num_fluid == 1 && variable_names[0].size() == 1)
	{
		tecplot_file << "VARLOCATION =([" << 1 + grid_data.DIM	<<"]=CELLCENTERED)" << endl;
	}
	else
	{
		int NV	= 0;
		for (int f = 0; f <= num_fluid-1; f++)	NV += variable_names[f].size();

		tecplot_file << "VARLOCATION =([" << grid_data.DIM+1 << "-" << NV+grid_data.DIM <<"]=CELLCENTERED)" << endl;
	}



	// WRITE DATA //
	// 1. NODE DATA
	// Write temporary files
	int node_ID_ptr;
	for (j = 0; j <= grid_data.DIM-1; j++)
	{
		for (i = 1; i <= grid_data.NNM; i++)
		{
			node_ID_ptr	= grid_data.nodes.whereis[i];
			tecplot_file << grid_data.nodes.data_ptr[node_ID_ptr]->x[j] << endl;
		}

		tecplot_file << endl;
	}


	// 2. FLOW DATA
	int cell_ID_ptr;
	for (int f = 0; f <= num_fluid-1; f++)
	{
		NV	= variable_names[f].size();
		for (j = 0; j <= NV-1; j++)
		{
			for (i = 1; i <= grid_data.NCM; i++)
			{
				cell_ID_ptr	= grid_data.cells.whereis[i];
				tecplot_file << Vc[f][cell_ID_ptr][j] <<endl;
			}
		}
	}


	// 3. WRITE CONNECTIVITY
	if (grid_data.DIM == 2)
	{
		for (i = 1; i <= grid_data.NCM; i++)
		{
			cell_ID_ptr	= grid_data.cells.whereis[i];

			switch (grid_data.cells.data_ptr[cell_ID_ptr]->type)
			{
			case C_QUADRILATERAL:
				for (j = 0; j <= 3; j++)
				{
					tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[j] << " ";
				}
				break;
			case C_TRIANGLE:
				for (j = 0; j <= 2; j++)
				{
					tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[j] << " ";
				}
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2];
				break;
			}

			tecplot_file << endl;
		}
	}
	else
	{
		for (i = 1; i <= grid_data.NCM; i++)
		{
			cell_ID_ptr	= grid_data.cells.whereis[i];
			switch (grid_data.cells.data_ptr[cell_ID_ptr]->type)
			{

			case C_TETRAHEDRON:
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[0] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[1] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3];
				tecplot_file << endl;
				break;

			case C_HEXAHEDRON:
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[0] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[1] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[5] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[6] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[7];
				tecplot_file << endl;
				break;

			case C_PYRAMID:
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[0] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[1] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4];
				tecplot_file << endl;
				break;

			case C_PRISM:
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[0] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[1] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[2] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[3] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[4] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[5] <<" ";
				tecplot_file << grid_data.cells.data_ptr[cell_ID_ptr]->node[6];
				tecplot_file << endl;
				break;
			}
		}
	} // END writing connectivity

	tecplot_file.close();
}







void OP2A_data_print_restart(int P, GRID_CLASS &grid_data, vector< vector < vector<double> > > &Vc, string file_name, int num_fluid, int iter, int type)
{
	string			filetype;
	ofstream		restart_file;

	// OPEN FILE TO WRITE
	restart_file.open(file_name.c_str());


	restart_file << "[ITERATION]: " << iter << endl;
	restart_file << "[NUMBER OF FLUID]:" << num_fluid << endl;
	restart_file << "[VARIABLES TYPE]:" << type << endl;
	restart_file << "[NCM]:" << grid_data.NCM << endl;

	for (int f = 0; f <= num_fluid-1; f++)
	{
		int NV	= Vc[f][0].size();
		for (int c = 1; c <= grid_data.NCM; c++)
		{
			int cell_ID_ptr	= grid_data.cells.whereis[c];
			for (int i = 0; i <=  NV-1; i++)
			{
				restart_file << Vc[f][cell_ID_ptr][i] << "  ";
			}

			restart_file << endl;
		}
	}


	restart_file.close();
}




