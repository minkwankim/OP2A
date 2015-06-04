/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 4, 2015
 *      			Author: Minkwan Kim
 *
 * GridReadFace.cpp
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


void read_mesh_face_fluent(const string& mesh_file_name,		// Mesh file name
							unsigned int	DIM,
							unsigned int	NFM,				// Grid information
							vector<Node>	&nodes,
							vector<Face>	&faces,
							vector<Cell>	&cells,
							vector<int>		&bc_zone)
{
	int index, index2, zone, first_index, last_index, bc_type, face_type, face_type2;
	int i, n;

	int type;
	int	nodes_num[4];
	int	cl, cr;


	string line;
	ifstream mesh_file(mesh_file_name.c_str());

	// Step 1: Open file to read
	if(!mesh_file) throw Common::ExceptionFileSystem (FromHere(), "Could not open file: " + mesh_file_name);

	// Step 2: Read Data
	n = 1;

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
				case FaceType::f_mixed:
					for (i = first_index; i <= last_index; i++)
					{
						face_type2 = -1;
						getline(mesh_file, line);
						sscanf(line.c_str(), "%d", &face_type2);

						switch (face_type2)
						{
						case FaceType::f_line: // LINE
							sscanf(line.c_str(), "%d %x %x %x %x", &type, &nodes_num[0], &nodes_num[1], &cl, &cr);

							faces[i].geo.allocate(DIM, FaceType::f_line);

							faces[i].geo.node_list[0] = &nodes[nodes_num[0]];
							faces[i].geo.node_list[1] = &nodes[nodes_num[1]];

							faces[i].geo.cl[0]	= &cells[cl];
							faces[i].geo.cr[0]	= &cells[cr];

							break;

						case FaceType::f_triangle:	// TRIANGLE
							sscanf(line.c_str(), "%d %x %x %x %x %x", &type, &nodes_num[0], &nodes_num[1], &nodes_num[2], &cl, &cr);

							faces[i].geo.allocate(DIM, FaceType::f_triangle);

							faces[i].geo.node_list[0] = &nodes[nodes_num[0]];
							faces[i].geo.node_list[1] = &nodes[nodes_num[1]];
							faces[i].geo.node_list[2] = &nodes[nodes_num[2]];

							faces[i].geo.cl[0]	= &cells[cl];
							faces[i].geo.cr[0]	= &cells[cr];
							break;

						case FaceType::f_quadrilateral: // QUADRILATERAL
							sscanf(line.c_str(), "%d %x %x %x %x %x %x", &type, &nodes_num[0], &nodes_num[1], &nodes_num[2], &nodes_num[3], &cl, &cr);

							faces[i].geo.allocate(DIM, FaceType::f_quadrilateral);

							faces[i].geo.node_list[0] = &nodes[nodes_num[0]];
							faces[i].geo.node_list[1] = &nodes[nodes_num[1]];
							faces[i].geo.node_list[2] = &nodes[nodes_num[2]];
							faces[i].geo.node_list[3] = &nodes[nodes_num[3]];

							faces[i].geo.cl[0]	= &cells[cl];
							faces[i].geo.cr[0]	= &cells[cr];

							break;

						default:
							throw ExceptionNotSupportedType (FromHere(), "It is not supported Face-type");
							break;
						}

						faces[i].geo.BC = bc_type;
						n++;
					}
					break;


				case FaceType::f_line: // LINE
					for (i = first_index; i <= last_index; i++)
					{
						getline(mesh_file, line);
						sscanf(line.c_str(), "%x %x %x %x", &nodes_num[0], &nodes_num[1], &cl, &cr);

						faces[i].geo.allocate(DIM, FaceType::f_line);

						faces[i].geo.node_list[0] = &nodes[nodes_num[0]];
						faces[i].geo.node_list[1] = &nodes[nodes_num[1]];

						faces[i].geo.cl[0]	= &cells[cl];
						faces[i].geo.cr[0]	= &cells[cr];

						faces[i].geo.BC 	= bc_type;
						n++;
					}
					break;

				case FaceType::f_triangle:	// TRIANGLE
					for (i = first_index; i <= last_index; i++)
					{
						getline(mesh_file, line);
						sscanf(line.c_str(), "%x %x %x %x %x", &nodes_num[0], &nodes_num[1], &nodes_num[2], &cl, &cr);

						faces[i].geo.allocate(DIM, FaceType::f_triangle);

						faces[i].geo.node_list[0] = &nodes[nodes_num[0]];
						faces[i].geo.node_list[1] = &nodes[nodes_num[1]];
						faces[i].geo.node_list[2] = &nodes[nodes_num[2]];


						faces[i].geo.cl[0]	= &cells[cl];
						faces[i].geo.cr[0]	= &cells[cr];

						faces[i].geo.BC = bc_type;
						n++;
					}
					break;

				case FaceType::f_quadrilateral: // QUADRILATERAL
					for (i = first_index; i <= last_index; i++)
					{
						getline(mesh_file, line);
						sscanf(line.c_str(), "%x %x %x %x %x %x", &nodes_num[0], &nodes_num[1], &nodes_num[2], &nodes_num[3], &cl, &cr);

						faces[i].geo.allocate(DIM, FaceType::f_triangle);

						faces[i].geo.node_list[0] = &nodes[nodes_num[0]];
						faces[i].geo.node_list[1] = &nodes[nodes_num[1]];
						faces[i].geo.node_list[2] = &nodes[nodes_num[2]];
						faces[i].geo.node_list[3] = &nodes[nodes_num[3]];

						faces[i].geo.cl[0]	= &cells[cl];
						faces[i].geo.cr[0]	= &cells[cr];

						faces[i].geo.BC = bc_type;
						n++;
					}
					break;

				default:
					throw ExceptionNotSupportedType (FromHere(), "It is not supported Face-type");
					break;
				}
			}// END READ FACE DATA
		}// END READ MESH FILE
	}

	//Step 2.2 Check Error
	if ((n-1) != NFM)
	{
		throw ExceptionGridDataMismatch (FromHere(), "PROBLEM IN Face DATA. TOTAL NUMBER OF Face DATA DOES NOT MATHCH WITH MESH INFOMATION DATA");
	}

}









}
}
