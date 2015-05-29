/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Feb 4, 2015
 *      			Author: Minkwan Kim
 *
 * problem_viscos.cpp
 * 			-  
 *  
 */






#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../include/OP2A_problem_viscous.hpp"
#include "../include/OP2A_setup_constants.hpp"


PROBLEM_VISCOUS::PROBLEM_VISCOUS()
{
	model_viscosity		= 0;
	model_conductivity	= 0;
	model_mixing_rule	= 1;

	viscous_relaxation_na	= 100;
	viscous_relaxation_nb	= 500;
	viscous_relaxation_min	= 0.01;
	viscous_relaxation_max	= 1.0;

	adiabatic		= false;
	catalytic		= false;
	radiative		= false;
	use_emissivity 	= false;

	Le				= 1.0;
}


PROBLEM_VISCOUS::~PROBLEM_VISCOUS()
{


}




void PROBLEM_VISCOUS::read_VISCOUS(string filename)
{
	int error_code;
	string line;
	string error_message;
	ifstream prob_file;


	int flag_adiabatic	= 0;
	int flag_catalytic	= 0;
	int	flag_radiative	= 0;
	int flag_emissivity	= 0;

	prob_file.open(filename.c_str());


	// READ PROBLEM SETUP DATA FILE
	if (prob_file.is_open())
	{
		while (! prob_file.eof())
		{
			getline(prob_file, line);
			data_read::read_line(line);

			// NUMERICAL METHOD SETTING - CFD
			data_read::get_data_from_string<int>(line,		"VISCOSITY MODEL: ", 0, model_viscosity);
			data_read::get_data_from_string<int>(line,		"THERMAL CONDUCTIVITY MODEL: ", 0, model_conductivity);
			data_read::get_data_from_string<int>(line,		"ADIABATIC WALL: ", 0, flag_adiabatic);
			data_read::get_data_from_string<int>(line,		"CATALYTIC WALL: ", 0, flag_catalytic);
			data_read::get_data_from_string<int>(line,		"RADIATIVE WALL: ", 0, flag_radiative);
			data_read::get_data_from_string<int>(line,		"USE EMISSIVITY: ", 0, flag_emissivity);

			data_read::get_data_from_string<int>(line,			"VISCOUS RELAXATION FUNCTION START: ", 	100, viscous_relaxation_na);
			data_read::get_data_from_string<int>(line,			"VISCOUS RELAXATION FUNCTION END: ", 	500, viscous_relaxation_nb);
			data_read::get_data_from_string<double>(line,		"VISCOUS RELAXATION FUNCTION MIN: ",	0.01, viscous_relaxation_min);
			data_read::get_data_from_string<double>(line,		"VISCOUS RELAXATION FUNCTION MAX: ",	1.00, viscous_relaxation_max);
		}
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "OP2A_problem_VISCOUS";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A PROBLEM SETUP DATA FILE. PLEASE CHECK THE FILE";
		error.print_message();
	}


	if (flag_adiabatic == 1)	adiabatic = true;
	else						adiabatic = false;

	if (flag_catalytic == 1)	catalytic = true;
	else						catalytic = false;

	if (flag_radiative == 1)	radiative = true;
	else						radiative = false;

	if (flag_emissivity == 1)	use_emissivity = true;
	else						use_emissivity = false;


	prob_file.close();
}
