/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: May 26, 2014
 *      			Author: Minkwan Kim
 *
 * grid.cpp
 * 			-  
 *  
 */


#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "../include/constants_grid.hpp"
#include "../include/node_class.hpp"
#include "../include/face_class.hpp"
#include "../include/cell_class.hpp"

#include "../include/read_grid.hpp"

#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../../MATRIX/include/vector_MK.hpp"

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
					 	int type)
{
	switch (type)
	{
	case FLUENT:
		read_mesh_info_fluent(mesh_file_name, DIM, NNM, NFM , NCM, NGM, bc_zone);
		break;
	default:
		Error_message_type	error;
		error.location_primary_name	= "read_grid.cpp";
		error.message = "IT IS NOT SUPPORTED MESH FILE TYPE. PLEASE CHECK MESH FILE!!!";
		error.print_message();
		break;
	}
}



/*
 * 1. Fluent format
 */
void read_mesh_info_fluent(string mesh_file_name,		// Mesh file name
							unsigned int& DIM,			// Dimenstion
							unsigned int& NNM,			// Total number of nodes in mesh
							unsigned int& NFM,			// Total number of faces in mesh
							unsigned int& NCM,			// Total number of cells in mesh
							unsigned int& NGM,			// Total number of ghost-cells in mesh
							vector<int>& bc_zone)		// BC Information
{
	int i;
	int index, temp, zone, bc_type;
	int zone_temp;
	int index_s, index_e;
	int nd, nfm, ncm, nnm, ngm;
	int nfm_zone, num_zone;
	int num_mesh_zone[MAX_ZONE_READ_MESH];

	char zone_type[20];
	string line;
	ifstream mesh_file;

	for (i = 0; i <= MAX_ZONE_READ_MESH-1; i++) num_mesh_zone[i]	= 0;

	// Step 1: Open file to read
	mesh_file.open(mesh_file_name.c_str());
	num_zone = 0;
	ngm = 0;


	//////////////////
	// Read file	//
	//////////////////
	if (mesh_file.is_open())
	{
		// Step 2: FIRST READ
		while (! mesh_file.eof())
		{
			index = -1;
			getline(mesh_file, line);
			sscanf(line.c_str(), "(%d ", &index);

			switch(index)
			{
			case FLU_INDEX_DIMENSIONS:
				nd = 0;
				sscanf(line.c_str(), "(%d %d)", &temp, &nd);
				DIM = nd;
				break;

			case FLU_INDEX_NODE:
				nnm = 0;
				index_s = 0;
				index_e = 0;
				sscanf(line.c_str(),"(%d (%x %x %x", &temp, &zone, &index_s, &index_e); // GETTING INDEX AND ZONE
				if (zone == 0) // READ GENERAL INFORMATION OF NODES
				{
					NNM = index_e - index_s + 1;
				}
				break;

			case FLU_INDEX_CELL:
				ncm = 0;
				index_s = 0;
				index_e = 0;
				sscanf(line.c_str(),"(%d (%x %x %x", &temp, &zone, &index_s, &index_e); // GETTING INDEX AND ZONE
				if (zone == 0) // READ GENERAL INFORMATION OF NODES
				{
					NCM = index_e - index_s + 1;
				}
				break;

			case FLU_INDEX_FACE:
				nfm = 0;
				index_s = 0;
				index_e = 0;
				sscanf(line.c_str(),"(%d (%x %x %x", &temp, &zone, &index_s, &index_e); // GETTING INDEX AND ZONE

				if (zone == 0) // READ GENERAL INFORMATION OF NODES
				{
					NFM = index_e - index_s + 1;
				}
				else
				{
					sscanf(line.c_str(),"(%d (%x %x %x %d", &temp, &zone_temp, &index_s, &index_e, &bc_type); // GETTING INDEX AND ZONE
					nfm_zone 			= index_e - index_s + 1;
					num_mesh_zone[zone] = nfm_zone;
					num_zone++;
				}
				break;

			case FLU_INDEX_BC:
				sscanf(line.c_str(),"(%d (%d %s", &index, &zone, zone_type); // GETTING INDEX AND ZONE

				if (zone > MAX_ZONE_READ_MESH)
					cerr << "[ERROR][READ_MESH]: PROBLEM READING MESH FILE - TOO MANY ZONES!! PLEASE ADJUST MAXIMUM LIMIT OF ZONES." << endl;

				// GENERAL TYPE OF BOUNDARY CONDITIONS
				if (strcmp(zone_type,"interior") == 0) 				bc_zone[zone] = BC_INTERIOR;
				if (strcmp(zone_type,"wall") == 0) 					bc_zone[zone] = BC_WALL;
				if (strcmp(zone_type,"pressure-inlet") == 0)		bc_zone[zone] = BC_INLET;
				if (strcmp(zone_type,"pressure-outlet") == 0)		bc_zone[zone] = BC_OUTLET;
				if (strcmp(zone_type,"pressure-far-field") == 0)	bc_zone[zone] = BC_FREESTREAM;
				if (strcmp(zone_type,"symmetry") == 0)				bc_zone[zone] = BC_SYMMETRY;
				if (strcmp(zone_type,"periodic-shadow") == 0)		bc_zone[zone] = BC_PERIODIC_SHADOW;
				if (strcmp(zone_type,"velocity-inlet") == 0)		bc_zone[zone] = BC_VELOCITY_INLET;
				if (strcmp(zone_type,"periodic") == 0)				bc_zone[zone] = BC_PERIODIC;
				if (strcmp(zone_type,"fan") == 0)					bc_zone[zone] = BC_FAN;
				if (strcmp(zone_type,"mass-flow-inlet") == 0)		bc_zone[zone] = BC_MASS_FLOW_INLET;
				if (strcmp(zone_type,"interface") == 0)				bc_zone[zone] = BC_INTERFACE;
				if (strcmp(zone_type,"parent") == 0)				bc_zone[zone] = BC_PARENT;
				if (strcmp(zone_type,"outflow") == 0)				bc_zone[zone] = BC_OUTFLOW;
				if (strcmp(zone_type,"axis") == 0)					bc_zone[zone] = BC_AXIS;
				if (strcmp(zone_type,"anode") == 0) 				bc_zone[zone] = BC_ANODE;
				if (strcmp(zone_type,"cathode") == 0) 				bc_zone[zone] = BC_CATHODE;
				if (strcmp(zone_type,"dielectric-wall") == 0) 		bc_zone[zone] = BC_DIELECTRIC_WALL;
				break;
			}
		}


		// Step 3. CHECK BOUNDARY CONDITION
		for (i = 0; i <= MAX_ZONE_READ_MESH-1; i++)
		{
			nfm_zone		= num_mesh_zone[i];
			if (bc_zone[i] != BC_INTERIOR)
			{
				ngm += nfm_zone;
			}
		}
		NGM = ngm;
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name	= "read_grid.cpp";
		error.message = "Cannot find grid file. Please check your grid file!";
		error.print_message();
	}

	mesh_file.close();
}







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
						int TYPE)
{
	switch (TYPE)
	{
	case FLUENT:
		read_mesh_node_fluent(mesh_file_name, DIM, NNM, nodes, whereis);
		break;

	default:
		Error_message_type	error;
		error.location_primary_name	= "read_grid.cpp";
		error.message = "THIS IS NOT SUPPORTED MESH FILE TYPE. PLEASE CHECK MESH FILE!!!";
		error.print_message();
		break;
	}
}


void read_mesh_node_fluent(string 		mesh_file_name,		// Mesh file name
							unsigned int	DIM,
							unsigned int	NNM,				// Grid information
							vector<NODE_CLASS *>	&nodes,
							vector<int>				&whereis)
{
	int i, n;
	int index, index2, zone, first_index, last_index, type, nd;
	double x, y, z;

	string line;
	ifstream mesh_file;

	whereis.resize(NNM+1);

	// Step 1: Open file to read
	mesh_file.open(mesh_file_name.c_str());
	n = 0;


	// Step 2: Read Data
	if (mesh_file.is_open())
	{
		while (! mesh_file.eof())
		{
			// STEP 1: GET INDEX
			index = -1;
			getline(mesh_file, line);
			sscanf(line.c_str(), "(%d ", &index);

			if (index == FLU_INDEX_NODE)
			{
				index2		= -1;
				zone 		= -1;
				first_index	= -1;
				last_index	= -1;
				type		= -1;
				nd			= -1;

				sscanf(line.c_str(),"(%d (%x %x %x %d %d)", &index2, &zone, &first_index, &last_index, &type, &nd); // GETTING INDEX AND ZONE
				if (zone > 0)
				{
					if (nd != DIM)
					{
						Error_message_type	error;
						error.location_primary_name	= "read_grid.cpp";
						error.message = "PROBLEM IN DIMENSION OF NODE DATA. Please the dimension of grid file";
						error.print_message();
					}
					else
					{
						switch (nd)
						{
						case 2:
							for (i = first_index; i <= last_index; i++)
							{
								x = 0.0;
								y = 0.0;

								getline(mesh_file, line);
								sscanf(line.c_str(), "%lf %lf", &x, &y);

								if (fabs(x) <= MATH_ZERO_MESH)	x = 0.0;
								if (fabs(y) <= MATH_ZERO_MESH)	y = 0.0;


								nodes.push_back(new NODE_CLASS());
								nodes[n]->ID	= i;
								nodes[n]->x[0]	= x;
								nodes[n]->x[1]	= y;

								whereis[i]	= n;
								n++;
							}
							break;
						case 3:
							for (i = first_index; i <= last_index; i++)
							{
								x = 0.0;
								y = 0.0;
								z = 0.0;

								getline(mesh_file, line);
								sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z);

								if (fabs(x) <= MATH_ZERO_MESH)	x = 0.0;
								if (fabs(y) <= MATH_ZERO_MESH)	y = 0.0;
								if (fabs(z) <= MATH_ZERO_MESH)	z = 0.0;

								nodes.push_back(new NODE_CLASS());
								nodes[n]->ID	= i;
								nodes[n]->x[0]	= x;
								nodes[n]->x[1]	= y;
								nodes[n]->x[2]	= z;

								whereis[i]	= n;

								n++;
							}
							break;
						}
					}
				}// END READ NODE DATA
			}

		}// END READ MESH FILE

		// Step 3: Check error
		if (n != NNM)
		{
			Error_message_type	error;
			error.location_primary_name	= "read_grid.cpp";
			error.message = "PROBLEM IN NODE DATA. TOTAL NUMBER OF NODE DATA DOES NOT MATHCH WITH MESH INFOMATION DATA";
			error.print_message();
		}
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name	= "read_grid.cpp";
		error.message = "CANNOT OPEN MESH FILE. PLEASE CHECK THE MESH FILE!!!";
		error.print_message();
	}
}







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
						int TYPE)
{
	switch (TYPE)
	{
	case FLUENT:
		read_mesh_face_fluent(mesh_file_name, DIM, NFM, faces, whereis, bc_zone);
		break;

	default:
		program_error_type1("THIS IS NOT SUPPORTED MESH FILE TYPE. PLEASE CHECK MESH FILE!!!");
		break;
	}
}



void read_mesh_face_fluent(string 		mesh_file_name,		// Mesh file name
							unsigned int	DIM,
							unsigned int	NFM,				// Grid information
							vector<FACE_CLASS *>	&faces,
							vector<int>				&whereis,
							vector<int>				&bc_zone)
{
	int index, index2, zone, first_index, last_index, bc_type, face_type, face_type2;
	int i, n;

	int type;
	int	nodes[4];
	int	cl, cr;


	string line;
	ifstream mesh_file;
	whereis.resize(NFM+1);

	// Step 1: Open file to read
	mesh_file.open(mesh_file_name.c_str());
	n = 0;


	// Step 2: Read Data
	if (mesh_file.is_open())
	{
		while (! mesh_file.eof())
		{
		// STEP 2.1: GET INDEX
			index = -1;
			getline(mesh_file, line);
			sscanf(line.c_str(), "(%d ", &index);

			if (index == FLU_INDEX_FACE)
			{
				index2			= -1;
				zone 			= -1;
				first_index 	= -1;
				last_index 		= -1;
				bc_type 		= -1;
				face_type 		= -1;

				sscanf(line.c_str(),"(%d (%x %x %x %d %d)", &index2, &zone, &first_index, &last_index, &bc_type, &face_type); // GETTING INDEX AND ZONE

				if (zone > 0)
				{
					// Assign BC type from data
					bc_type	= bc_zone[zone];

					switch (face_type)
					{
					case F_MIXED:
						for (i = first_index; i <= last_index; i++)
						{
							faces.push_back(new FACE_CLASS());

							face_type2 = -1;
							getline(mesh_file, line);
							sscanf(line.c_str(), "%d", &face_type2);

							switch (face_type2)
							{
							case F_LINE: // LINE
								sscanf(line.c_str(), "%d %x %x %x %x", &type, &nodes[0], &nodes[1], &cl, &cr);

								faces[n]->ID	= i;
								faces[n]->type	= type;
								faces[n]->NN	= 2;
								faces[n]->node.resize(2);
								faces[n]->node[0]	= nodes[0];
								faces[n]->node[1]	= nodes[1];
								faces[n]->cl[0]	= cl;
								faces[n]->cr[0]	= cr;

								break;

							case F_TRIANGLE:	// TRIANGLE
								sscanf(line.c_str(), "%d %x %x %x %x %x", &type, &nodes[0], &nodes[1], &nodes[2], &cl, &cr);

								faces[n]->ID	= i;
								faces[n]->type	= type;
								faces[n]->NN	= 3;
								faces[n]->node.resize(3);
								faces[n]->node[0]	= nodes[0];
								faces[n]->node[1]	= nodes[1];
								faces[n]->node[2]	= nodes[2];
								faces[n]->cl[0]	= cl;
								faces[n]->cr[0]	= cr;

								break;

							case F_QUADRILATERAL: // QUADRILATERAL
								sscanf(line.c_str(), "%d %x %x %x %x %x %x", &type, &nodes[0], &nodes[1], &nodes[2], &nodes[3], &cl, &cr);

								faces[n]->ID	= i;
								faces[n]->type	= type;
								faces[n]->NN	= 4;
								faces[n]->node.resize(4);
								faces[n]->node[0]	= nodes[0];
								faces[n]->node[1]	= nodes[1];
								faces[n]->node[2]	= nodes[2];
								faces[n]->node[3]	= nodes[3];
								faces[n]->cl[0]	= cl;
								faces[n]->cr[0]	= cr;

								break;

							default:
								cout <<  "Face-type number" << face_type2 << " IS NOT SUPPORTED FACE TYPE" << endl;
								program_error_type1("PLEASE CHECK THE FACE TYPE IN A MESH FILE!");
								break;
							}

							faces[n]->BC = bc_type;
							whereis[i]	= n;

							n++;
						}
						break;

					case F_LINE: // LINE
						for (i = first_index; i <= last_index; i++)
						{
							faces.push_back(new FACE_CLASS());

							getline(mesh_file, line);
							sscanf(line.c_str(), "%x %x %x %x", &nodes[0], &nodes[1], &cl, &cr);

							faces[n]->ID	= i;
							faces[n]->type	= face_type;
							faces[n]->NN	= 2;
							faces[n]->node.resize(2);
							faces[n]->node[0]	= nodes[0];
							faces[n]->node[1]	= nodes[1];
							faces[n]->cl[0]	= cl;
							faces[n]->cr[0]	= cr;
							faces[n]->BC 	= bc_type;

							whereis[i]	= n;
							n++;
						}
						break;

					case F_TRIANGLE:	// TRIANGLE
						for (i = first_index; i <= last_index; i++)
						{
							faces.push_back(new FACE_CLASS());

							getline(mesh_file, line);
							sscanf(line.c_str(), "%x %x %x %x %x", &nodes[0], &nodes[1], &nodes[2], &cl, &cr);

							faces[n]->ID	= i;
							faces[n]->type	= face_type;
							faces[n]->NN	= 3;
							faces[n]->node.resize(3);
							faces[n]->node[0]	= nodes[0];
							faces[n]->node[1]	= nodes[1];
							faces[n]->node[2]	= nodes[2];
							faces[n]->cl[0]	= cl;
							faces[n]->cr[0]	= cr;
							faces[n]->BC = bc_type;

							whereis[i]	= n;
							n++;
						}
						break;

					case F_QUADRILATERAL: // QUADRILATERAL
						for (i = first_index; i <= last_index; i++)
						{
							faces.push_back(new FACE_CLASS());

							getline(mesh_file, line);
							sscanf(line.c_str(), "%x %x %x %x %x %x", &nodes[0], &nodes[1], &nodes[2], &nodes[2], &cl, &cr);

							faces[n]->ID	= i;
							faces[n]->type	= face_type;
							faces[n]->NN	= 4;
							faces[n]->node.resize(4);
							faces[n]->node[0]	= nodes[0];
							faces[n]->node[1]	= nodes[1];
							faces[n]->node[2]	= nodes[2];
							faces[n]->node[3]	= nodes[3];
							faces[n]->cl[0]	= cl;
							faces[n]->cr[0]	= cr;
							faces[n]->BC = bc_type;

							whereis[i]	= n;
							n++;
						}
						break;

					default:
						cout <<  "Face-type number" << face_type2 << " IS NOT SUPPORTED FACE TYPE" << endl;
						program_error_type1("PLEASE CHECK THE FACE TYPE IN A MESH FILE!");
						break;
					}
				}// END READ FACE DATA
			}// END READ MESH FILE
		}

		//Step 2.2 Check Error
		if (n != NFM)
		{
			cout << "PROBLEM IN FACE DATA. TOTAL NUMBER OF FACE DATA DOES NOT MATHCH WITH MESH INFOMATION DATA" << endl;
			program_error_type1("Please check grid file");
		}
	}
	else
	{
		program_error_type1("CANNOT OPEN MESH FILE. PLEASE CHECK THE MESH FILE!!!");
	}
}









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
						int TYPE)
{
	switch (TYPE)
	{
	case FLUENT:
		read_mesh_cell_fluent(mesh_file_name, DIM, NCM, cells, whereis);
		break;

	default:
		program_error_type1("T IS NOT SUPPORTED MESH FILE TYPE. PLEASE CHECK MESH FILE!!!");
		break;
	}
}


void read_mesh_cell_fluent(string mesh_file_name,		// Mesh file name
		unsigned int	DIM,
		unsigned int	NCM,				// Grid information
		vector<CELL_CLASS *>	&cells,
		vector<int>				&whereis)
{
	int i, n;
	int index, index2, zone, first_index, last_index, type, element_type;

	string line;
	ifstream mesh_file;
	whereis.resize(NCM+1);


	// Step 1: Open grid file
	mesh_file.open(mesh_file_name.c_str());
	n = 0;


	// Step 2: Read mesh file
	if (mesh_file.is_open())
	{
		while (! mesh_file.eof())
		{
			// STEP 2.1: GET INDEX
			index = -1;
			getline(mesh_file, line);
			sscanf(line.c_str(), "(%d ", &index);

			if (index == FLU_INDEX_CELL)
			{
				index2			= -1;
				zone 			= -1;
				first_index		= -1;
				last_index		= -1;
				type 			= -1;
				element_type	= -1;

				sscanf(line.c_str(),"(%d (%x %x %x %d", &index2, &zone, &first_index, &last_index, &type); // GETTING INDEX AND ZONE
				if (zone > 0)
				{
					sscanf(line.c_str(),"(%d (%x %x %x %d %d)", &index2, &zone, &first_index, &last_index, &type, &element_type); // GETTING INDEX AND ZONE

					if (type == 1)
					{
						if (element_type != C_MIXED)
						{
							for (i = first_index; i <= last_index; i++)
							{
								cells.push_back(new CELL_CLASS());

								cells[n]->ID	= i;
								cells[n]->BC	= type;
								cells[n]->type	= element_type;

								whereis[i]	= n;
								n++;
							}
						}
						else
						{
							for (i = first_index; i <= last_index; i++)
							{
								mesh_file >> element_type;
								cells.push_back(new CELL_CLASS());

								cells[n]->ID	= i;
								cells[n]->BC	= type;
								cells[n]->type	= element_type;

								whereis[i]	= n;
								n++;
							}
						}
					}
				}
			} // END READ CELL DATA
		} // END READ MESH FILE

		if (n != NCM)
		{
			program_error_type1("PROBLEM IN CELL DATA. TOTAL NUMBER OF CELL DATA DOES NOT MATHCH WITH MESH INFOMATION DATA");
		}
	}
	else
	{
		program_error_type1("CANNOT OPEN MESH FILE. PLEASE CHECK THE MESH FILE!!!");
	}
}
