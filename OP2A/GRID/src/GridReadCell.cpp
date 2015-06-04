/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 4, 2015
 *      			Author: Minkwan Kim
 *
 * GridReadCell.cpp
 * 			-  
 *  
 */




#include "Common/include/Exception_FileSystem.hpp"
#include "GRID/include/GridRead.hpp"
#include "GRID/include/Exception_DimensionMismatch.hpp"
#include "GRID/include/Exception_GridDataMismatch.hpp"
#include "GRID/include/Exception_NotSupportedType.hpp"

namespace OP2A{
namespace GRID{


void read_mesh_cell_fluent(const string& mesh_file_name,		// Mesh file name
							unsigned int	DIM,
							unsigned int	NCM,				// Grid information
							vector<Cell>	&cells,
							vector<int>		&bc_zone)
{

	int i, n;
	int index, index2, zone, first_index, last_index, type, element_type;

	string line;
	ifstream mesh_file(mesh_file_name.c_str());

	// Step 1: Open file to read
	if(!mesh_file) throw Common::ExceptionFileSystem (FromHere(), "Could not open file: " + mesh_file_name);


	// Step 2: Read mesh file
	n = 1;

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

				if (type != 0)
				{
					if (element_type != CellType::c_mixed)
					{
						for (i = first_index; i <= last_index; i++)
						{
							cells[i].geo.allocate(DIM, static_cast<CellType>(element_type));
							cells[i].geo.BC	= bc_zone[zone];
							n++;
						}
					}
					else
					{
						for (i = first_index; i <= last_index; i++)
						{
							mesh_file >> element_type;
							cells[i].geo.allocate(DIM, static_cast<CellType>(element_type));
							cells[i].geo.BC	= bc_zone[zone];
							n++;
						}
					}
				}
			}
		} // END READ CELL DATA
	} // END READ MESH FILE

	if ((n-1) != NCM)
	{
		throw ExceptionGridDataMismatch (FromHere(), "PROBLEM IN Cell DATA. TOTAL NUMBER OF Cell DATA DOES NOT MATHCH WITH MESH INFOMATION DATA");
	}

}








}
}
