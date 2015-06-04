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
}




}
}
