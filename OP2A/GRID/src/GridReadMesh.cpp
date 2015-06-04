/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 4, 2015
 *      			Author: Minkwan Kim
 *
 * GridReadMesh.cpp
 * 			-  
 *  
 */



#include "Common/include/Exception_FileSystem.hpp"

#include "GRID/include/Grid.hpp"
#include "GRID/include/GridRead.hpp"
#include "GRID/include/Exception_DimensionMismatch.hpp"
#include "GRID/include/Exception_GridDataMismatch.hpp"
#include "GRID/include/Exception_NotSupportedType.hpp"

namespace OP2A{
namespace GRID{

void Grid::readMeahData(const string& mesh_file_name, GridDataType type)
{
	 readMeahDataInfo(mesh_file_name, type);
	 readMeahDataNode(mesh_file_name, type);
	 readMeahDataFace(mesh_file_name, type);
	 readMeahDataCell(mesh_file_name, type);
}




}
}
