/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 29, 2015
 *      			Author: Minkwan Kim
 *
 * Grid.hpp
 * 			-  
 *  
 */
#ifndef GRID_HPP_
#define GRID_HPP_


#include "Node.hpp"
#include "Face.hpp"
#include "Cell.hpp"
#include "GridConstant.hpp"

#include "Common/include/Map1D.hpp"
#include "Common/include/Map2D.hpp"
#include "Common/include/Map3D.hpp"


namespace OP2A{
namespace GRID{

class Grid
{
public:
	unsigned int ND;
	unsigned int NNM;
	unsigned int NFM;
	unsigned int NCM;
	unsigned int NGM;

	vector<Node>	nodes;
	vector<Face>	faces;
	vector<Cell>	cells;
	vector<Cell>	cells_ghost;

public:	// Constructor and desctructor
	Grid();
	~Grid();

private:
	vector<int>	bc_zone;
	bool is_read_basic_info;

public:	// Functions
	void readMeshData(const string& mesh_file_name, GridDataType type);
	void processingGridData(const double mesh_factor, bool is_axisymmetric);



private:
	void readMeshDataInfo(const string& mesh_file_name, GridDataType type);
	void readMeshDataNode(const string& mesh_file_name, GridDataType type);
	void readMeshDataFace(const string& mesh_file_name, GridDataType type);
	void readMeshDataCell(const string& mesh_file_name, GridDataType type);

	void processingNodeData(const double mesh_factor, bool is_axisymmetric);
	void processingFaceData();
	void processingCellData();
	void processingGhostData();


};







}
}



#endif /* GRID_HPP_ */
