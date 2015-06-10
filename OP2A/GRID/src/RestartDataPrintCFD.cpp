/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 10, 2015
 *      			Author: Minkwan Kim
 *
 * RestartDataPrintCFD.cpp
 * 			-  
 *  
 */



#include <fstream>
#include "GRID/include/PrintResult.hpp"


namespace OP2A{
namespace GRID{


void RestartDataPrintCFD(int P, unsigned int iter, Grid& grid, const string& file_name, string& data_used_to_restart)
{
	ofstream restart_file(file_name.c_str());

	int data_index = grid.cells[1].data.dataMap.find(data_used_to_restart);
	int NV	= grid.cells[1].data.data[data_index].numData;


	restart_file << "[ITERATION]: " << iter << endl;
	restart_file << "[NUMBER OF VARIABLES]:" << NV << endl;
	restart_file << "[NCM]:" << grid.NCM << endl;

	for (int c = 1; c <= grid.NCM; c++)
	{
		restart_file << grid.cells[c].geo.ID << "  ";

		for (int i = 0; i <=  NV-1; i++)
		{
			restart_file << grid.cells[c].data.data[data_index].data[i] << "  ";
		}

		restart_file << endl;
	}

	restart_file.close();
}





}
}
