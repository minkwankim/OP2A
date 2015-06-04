/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 4, 2015
 *      			Author: Minkwan Kim
 *
 * Grid.cpp
 * 			-  
 *  
 */



#include <vector>
using namespace std;

#include "../include/Grid.hpp"
#include "../include/GridRead.hpp"


namespace OP2A{
namespace GRID{

Grid::Grid():ND(2), NNM(0), NFM(0), NCM(0), NGM(0), is_read_basic_info(false)
{

}

Grid::~Grid()
{

}

/*
 * Internal functions
 */
void Grid::readMeahDataInfo(const string& mesh_file_name, GridDataType type)
{
	switch (type)
	{
	case GridDataType::FLUENT:
		read_mesh_info_fluent(mesh_file_name, ND, NNM, NFM,	NCM, NGM, bc_zone);
	break;
	}

	is_read_basic_info	= true;

	nodes.resize(NNM+1);
	faces.resize(NFM+1);
	cells.resize(NCM+1);
	cells_ghost.resize(NGM+1);
}


void Grid::readMeahDataNode(const string& mesh_file_name, GridDataType type)
{
	if (is_read_basic_info == false)	readMeahDataInfo(mesh_file_name, type);

	switch (type)
	{
	case GridDataType::FLUENT:
		read_mesh_node_fluent(mesh_file_name, ND, NNM, nodes);
	break;
	}
}


void Grid::readMeahDataFace(const string& mesh_file_name, GridDataType type)
{
	if (is_read_basic_info == false)	readMeahDataInfo(mesh_file_name, type);

	switch (type)
	{
	case GridDataType::FLUENT:
		read_mesh_face_fluent(mesh_file_name, ND, NFM, nodes, faces, cells, bc_zone);
	break;
	}
}

void Grid::readMeahDataCell(const string& mesh_file_name, GridDataType type)
{
	if (is_read_basic_info == false)	readMeahDataInfo(mesh_file_name, type);

	switch (type)
	{
	case GridDataType::FLUENT:
		read_mesh_cell_fluent(mesh_file_name, ND, NCM, cells, bc_zone);
	break;
	}
}









}
}
