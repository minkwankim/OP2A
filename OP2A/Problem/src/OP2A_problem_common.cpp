/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 13, 2015
 *      			Author: Minkwan Kim
 *
 * problem_common.cpp
 * 			-  
 *  
 */

#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../include/OP2A_problem_common.hpp"



/*
 * COMMON PROBLEM SETUP CLASS
 */
PROBLEM_COMMON::PROBLEM_COMMON()
{
	name	= "OPPA_default";						/* PROBLEM NAME				*/

	mesh_file_name	="OPPA_mesh.cas";
	mesh_file_type	= 0;
	grid_factor		= 1.0;

	output_file_name	= "OPPA_resut.plt";
	output_file_type	= 0;					/* OUTPUT FILE FORMAT				*/
	use_multi_file		= false;					/* Multi file						*/

	itv_save			= 100;						/* Interval for writing restart file				*/
	itv_result			= 100;						/* Interval for writing result file					*/

	n_current			= 0;						// Current simulation time step
	n_total				= 1e5;            			// Maximum number of time steps
	n_startup			= 10000;        			// Number of time steps until stationary
	convergence_criterion	= 1.0e-10;

	is_axisymmetric		= false;
	is_viscous			= false;
	DIM					= 2;
	NE					= 1;
	NCM					= 0;

	NP					= 1;						/* Number of CPUs	*/
	NT					= 4;
	P					= 0;						// Current CPU ID;
};

PROBLEM_COMMON::~PROBLEM_COMMON()
{

}




void PROBLEM_COMMON::read_COMMON(string file_name)
{
	int num;
	int error_code;
	int i_index, j_index;
	int NS_temp, ND_temp, NE_temp;
	int index_temp;
	int s_temp;
	int axis_flag;
	int viscous_flag;
	string line;
	string string_temp;
	string temp_text;
	string temp_text1;
	string temp_text2;
	string temp_text3;
	double temp1;
	double temp2;
	bool 	flag1;

	string error_message;
	ifstream prob_file;

	prob_file.open(file_name.c_str());


	// READ PROBLEM SETUP DATA FILE
	if (prob_file.is_open())
	{
		while (! prob_file.eof())
		{
			getline(prob_file, line);
			data_read::read_line(line);

			// BASIC INFOMATION
			string_temp = data_read::get_data_from_string_string(line, "TITLE OF PROBLEM:");
			if(string_temp.size() > 0)	name	= string_temp;


			string_temp = data_read::get_data_from_string_string(line, "MESH FILE NAME:");
			if(string_temp.size() > 0)	mesh_file_name	= string_temp;
			data_read::get_data_from_string<unsigned int>(line,	"MESH FILE TYPE:", 0, mesh_file_type);
			data_read::get_data_from_string<double>(line,	"GRID FACTOR:", 0, grid_factor);


			string_temp = data_read::get_data_from_string_string(line, "OUTPUT FILE NAME:");
			if(string_temp.size() > 0)	output_file_name	= string_temp;
			data_read::get_data_from_string<unsigned int>(line,	"OUTPUT FILE TYPE:", output_file_type, output_file_type);
			data_read::get_data_from_string<bool>(line,	"MULTI FILE OPTION:", false, use_multi_file);
			data_read::get_data_from_string<unsigned int>(line,		"INTERVAL FOR WRITING RESTART DATA:", 100, itv_save);
			data_read::get_data_from_string<unsigned int>(line,		"INTERVAL FOR WRITING RESULT DATA:", 100, itv_result);


			data_read::get_data_from_string<unsigned int>(line,		"MAXIMUM NUMBER OF TIME STEP:", 1e5, n_total);
			data_read::get_data_from_string<unsigned int>(line,		"NUMBER OF TIME STEP UNTIL STATIONARY:", 100, n_startup);
			data_read::get_data_from_string<double>(line,		"CONVERGENCE CRETERION:", 1.0e-10, convergence_criterion);


			data_read::get_data_from_string<int>(line,	"IS AXISYMMETRIC SIMULATION:", 0, axis_flag);
			data_read::get_data_from_string<int>(line,	"VISCOUS FLOW:", 0, viscous_flag);

		}

		if (axis_flag == 1)	is_axisymmetric	= true;
		else 				is_axisymmetric	= false;

		if (viscous_flag	== 1)	is_viscous	= true;
		else						is_viscous	= false;


		if (!name.empty() && name[name.size() - 1] == '\r')	name.erase(name.size() - 1);
		if (!mesh_file_name.empty())
		{
			if (mesh_file_name[mesh_file_name.size() - 1] == '\r' 	|| mesh_file_name[mesh_file_name.size() - 1] == '\t')
				mesh_file_name.erase(mesh_file_name.size() - 1);
		}

		if (!output_file_name.empty())
		{
			if (output_file_name[output_file_name.size() - 1] == '\r' 	|| output_file_name[output_file_name.size() - 1] == '\t')
				output_file_name.erase(output_file_name.size() - 1);
		}
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "OP2A_problem_common";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A PROBLEM SETUP DATA FILE. PLEASE CHECK THE FILE";
		error.print_message();
	}

	prob_file.close();
}
