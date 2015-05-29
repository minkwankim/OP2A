/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Feb 4, 2015
 *      			Author: Minkwan Kim
 *
 * problem_NONEQ.cpp
 * 			-  
 *  
 */



#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../include/OP2A_problem_NONEQ.hpp"
#include "../include/OP2A_setup_constants.hpp"


PROBLEM_NONEQ::PROBLEM_NONEQ()
{
	NER			= false;
	NEV			= false;
	NEE			= false;
	NONEQ_Chem	= false;

	ID_Tr	= 0;
	ID_Tv	= 0;
	ID_Te	= 0;
}

PROBLEM_NONEQ::~PROBLEM_NONEQ()
{

}


void PROBLEM_NONEQ::read_NONEQ(string filename)
{
	int error_code;
	string line;
	string error_message;
	string string_temp;
	ifstream prob_file;

	int ner, nev, nee, noneq_chem;

	ner	= 0;
	nev = 0;
	nee = 0;
	noneq_chem = 0;

	prob_file.open(filename.c_str());


	// READ PROBLEM SETUP DATA FILE
	if (prob_file.is_open())
	{
		while (! prob_file.eof())
		{
			getline(prob_file, line);
			data_read::read_line(line);

			// NUMERICAL METHOD SETTING - CFD
			data_read::get_data_from_string<int>(line,		"ROTATIONAL TEMPERATURE NONEQULIBRIUM:", 0, NER);
			data_read::get_data_from_string<int>(line,		"VIBRATIONAL TEMPERATURE NONEQULIBRIUM:", 0, NEV);
			data_read::get_data_from_string<int>(line,		"ELECTRON TEMPERATURE NONEQULIBRIUM:", 0, NEE);
			data_read::get_data_from_string<int>(line,		"CHEMICAL NONEQUILIBRIUM: ", 0, noneq_chem);

			string_temp = data_read::get_data_from_string_string(line, "CHEMICAL REACTION FILE:");
			if(string_temp.size() > 0)
			{
				reaction_file	= string_temp;
			}
		}

		if (noneq_chem == 1) NONEQ_Chem = true;
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "OP2A_problem_NONEQ";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A PROBLEM SETUP DATA FILE. PLEASE CHECK THE FILE";
		error.print_message();
	}

	prob_file.close();
}
