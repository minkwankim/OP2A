/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 9, 2015
 *      			Author: Minkwan Kim
 *
 * PrintResult.hpp
 * 			-  
 *  
 */
#ifndef PRINTRESULT_HPP_
#define PRINTRESULT_HPP_


#include "Grid.hpp"

namespace OP2A{
namespace GRID{

void ResultDataPrintTecplotCell(int P, Grid& grid, const string& title, const string& file_name, const string& data_to_print);
void ResultDataPrintTecplotNode(int P, Grid& grid, const string& title, const string& file_name, const string& data_to_print);
void ResultDataPrintTecplotNode2D(int P, Grid& grid, const string& title, const string& file_name, const string& data_to_print);
void ResultDataPrintTecplotNode3D(int P, Grid& grid, const string& title, const string& file_name, const string& data_to_print);

void RestartDataPrintCFD(int P, unsigned int iter, Grid& grid, const string& file_name, const string& data_used_to_restart);



}
}
#endif /* PRINTRESULT_HPP_ */
