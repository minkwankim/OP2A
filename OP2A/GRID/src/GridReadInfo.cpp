/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 4, 2015
 *      			Author: Minkwan Kim
 *
 * GridReadInfo.cpp
 * 			-  
 *  
 */


#include "Common/include/Exception_FileSystem.hpp"
#include "GRID/include/GridRead.hpp"


namespace OP2A{
namespace GRID{


void read_mesh_info_fluent(const string& mesh_file_name,	// Mesh file name
							unsigned int& DIM,				// Dimenstion
							unsigned int& NNM,				// Total number of nodes in mesh
							unsigned int& NFM,				// Total number of faces in mesh
							unsigned int& NCM,				// Total number of cells in mesh
							unsigned int& NGM,				// Total number of ghost-cells in mesh
							vector<int>& bc_zone)			// BC Information
{
	int i;
	int index, temp, zone, bc_type;
	int zone_temp;
	int index_s, index_e;
	int nd, nfm, ncm, nnm, ngm;
	int nfm_zone, num_zone;
	vector<int> num_mesh_zone(MAX_ZONE_READ_MESH, 0);
	bc_zone.resize(MAX_ZONE_READ_MESH, 0);

	char zone_type[20];
	string line;
	ifstream mesh_file(mesh_file_name.c_str());


	// Step 1: Open file to read
	if(!mesh_file) throw Common::ExceptionFileSystem (FromHere(), "Could not open file: " + mesh_file_name);

	num_zone = 0;
	ngm = 0;


	//////////////////
	// Read file	//
	//////////////////

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
			if (strcmp(zone_type,"interior") == 0) 				bc_zone[zone] = BCType::interior;
			if (strcmp(zone_type,"wall") == 0) 					bc_zone[zone] = BCType::wall;
			if (strcmp(zone_type,"pressure-inlet") == 0)		bc_zone[zone] = BCType::inlet;
			if (strcmp(zone_type,"pressure-outlet") == 0)		bc_zone[zone] = BCType::outlet;
			if (strcmp(zone_type,"pressure-far-field") == 0)	bc_zone[zone] = BCType::freestream;
			if (strcmp(zone_type,"symmetry") == 0)				bc_zone[zone] = BCType::symmetric;
			if (strcmp(zone_type,"axis") == 0)					bc_zone[zone] = BCType::axis;
			if (strcmp(zone_type,"anode") == 0) 				bc_zone[zone] = BCType::anode;
			if (strcmp(zone_type,"cathode") == 0) 				bc_zone[zone] = BCType::cathode;
			if (strcmp(zone_type,"dielectric-wall") == 0) 		bc_zone[zone] = BCType::dielectricwall;
			break;
		}
	}


	// Step 3. CHECK BOUNDARY CONDITION
	for (i = 0; i <= MAX_ZONE_READ_MESH-1; i++)
	{
		nfm_zone		= num_mesh_zone[i];
		if (bc_zone[i] != BCType::interior)
		{
			ngm += nfm_zone;
		}
	}
	NGM = ngm;

	mesh_file.close();
}

}
}

