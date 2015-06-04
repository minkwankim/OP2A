/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 4, 2015
 *      			Author: Minkwan Kim
 *
 * GridReadNode.cpp
 * 			-  
 *  
 */



#include "Common/include/Exception_FileSystem.hpp"
#include "GRID/include/GridRead.hpp"
#include "GRID/include/Exception_DimensionMismatch.hpp"
#include "GRID/include/Exception_GridDataMismatch.hpp"


namespace OP2A{
namespace GRID{


void read_mesh_node_fluent(const string& mesh_file_name,		// Mesh file name
							unsigned int	DIM,
							unsigned int	NNM,				// Grid information
							vector<Node>	&nodes)
{
	int i, n;
	int index, index2, zone, first_index, last_index, type, nd;
	double x, y, z;

	string line;
	ifstream mesh_file(mesh_file_name.c_str());

	// Step 1: Open file to read
	if(!mesh_file) throw Common::ExceptionFileSystem (FromHere(), "Could not open file: " + mesh_file_name);

	// Step 2: Read Data
	n = 0;
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
				if (nd != DIM)	throw ExceptionDimensionMismatch (FromHere(), "Problem in dimension of node data");

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

						nodes[n].geo.x[0]	= x;
						nodes[n].geo.x[1]	= y;
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

						nodes[n].geo.x[0]	= x;
						nodes[n].geo.x[1]	= y;
						nodes[n].geo.x[2]	= z;
						n++;
					}
					break;
				}
			}// END READ NODE DATA
		}
	}

	// Step 3: Check error
	if (n != NNM)
	{
		throw ExceptionGridDataMismatch (FromHere(), "PROBLEM IN NODE DATA. TOTAL NUMBER OF NODE DATA DOES NOT MATHCH WITH MESH INFOMATION DATA");
	}
}


}
}
