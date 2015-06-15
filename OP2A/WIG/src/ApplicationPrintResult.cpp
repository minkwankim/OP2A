/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 15, 2015
 *      			Author: Minkwan Kim
 *
 * ApplicationPrintResult.cpp
 * 			-  
 *  
 */


#include "Common/include/CodeLocation.hpp"
#include "Common/include/Exception_NPExceed.hpp"
#include "Common/include/StringOps.hpp"

#include "../include/OP2A_Application.hpp"


void ApplicationOP2A::print_result()
{
	string print_data_string = "Primitive variables";
	GRID::ResultDataPrintTecplotCell(P, grid, problem_setup.name, problem_setup.output_file_name, print_data_string);

	check_elapsed_time("Print Solution Data");
}

void ApplicationOP2A::print_restartCFD()
{
	string data_used_to_restart = "Primitive variables";
	string restart_file_name = problem_setup.name + ".rst";

	GRID::RestartDataPrintCFD(P, iter, grid, restart_file_name, data_used_to_restart);

	check_elapsed_time("Save Restart Data");
}
