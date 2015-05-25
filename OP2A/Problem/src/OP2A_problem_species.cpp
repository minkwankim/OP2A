/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 14, 2015
 *      			Author: Minkwan Kim
 *
 * problem_species.cpp
 * 			-  
 *  
 */


#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../include/OP2A_problem_species.hpp"
#include "../include/OP2A_setup_constants.hpp"


PROBLEM_SPECIES_BASIC::PROBLEM_SPECIES_BASIC()
{
	NS	= 0;
}

PROBLEM_SPECIES_BASIC::~PROBLEM_SPECIES_BASIC()
{

}

void PROBLEM_SPECIES_BASIC::read_SPECIES_BASIC(string file_name)
{
	int error_code;
	string string_temp;
	string line;
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

			data_read::get_data_from_string<int>(line,	"NUMBER OF SPECIES:", 0, NS);

			string_temp = data_read::get_data_from_string_string(line, "SPECIES FILE:");
			if(string_temp.size() > 0)
			{
				species_file	= string_temp;
			}
		}
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "OP2A_problem_SPECIES";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A PROBLEM SETUP DATA FILE. PLEASE CHECK THE FILE";
		error.print_message();
	}

	prob_file.close();
}
