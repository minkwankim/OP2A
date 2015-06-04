/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 4, 2015
 *      			Author: Minkwan Kim
 *
 * GridRead.hpp
 * 			-  
 *  
 */
#ifndef GRIDREAD_HPP_
#define GRIDREAD_HPP_


#include <vector>
#include <cmath>
#include <string.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "GridConstant.hpp"


using namespace std;


namespace OP2A{
namespace GRID{


#define MAX_ZONE_READ_MESH	1000

/*
 * 1. Fluent format
 */
void read_mesh_info_fluent(const string& mesh_file_name,		// Mesh file name
							unsigned int& DIM,					// Dimenstion
							unsigned int& NNM,					// Total number of nodes in mesh
							unsigned int& NFM,					// Total number of faces in mesh
							unsigned int& NCM,					// Total number of cells in mesh
							unsigned int& NGM,					// Total number of ghost-cells in mesh
							vector<int>& bc_zone);				// BC Information




}
}


#endif /* GRIDREAD_HPP_ */
