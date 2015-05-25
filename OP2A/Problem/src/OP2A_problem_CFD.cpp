/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 13, 2015
 *      			Author: Minkwan Kim
 *
 * problem_CFD.cpp
 * 			-  
 *  
 */


#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../include/OP2A_problem_CFD.hpp"
#include "../include/OP2A_setup_constants.hpp"


PROBLEM_CFD::PROBLEM_CFD()
{
	multi_fluid						= 1;
	SPATIAL_INTEGRATION_METHOD		= 0;			/* SPARTIAL INTEGRATION METHOD		*/
	TIME_INTEGRATION_METHOD			= 0;				/* TIME INTEGRATION METHOD			*/
	NUMERICAL_ORDER					= 1;					/* NUMERICAL ORDER					*/
	LIMITER							= 0;					/* LIMITER METHODS					*/
	USE_least_square_method			= 0;

	CFL_INCREASE_METHOD				= 0;					/* CFL NUMBER INCREASE METHOD (IT IS ONLY USED WHEN IMPLICIT METHODS ARE EMPLOYED)	*/
	iteration_before_1				= 1000;					/* NUMBER OF INTERATION BEFORE CFL NUMBER LARGER THAN 1 IS APPLIED					*/
	CFL_start						= 0.0001;				/* INITIAL CFL NUMBER (FOR EXPLICIT METHOD, IT IS AN APPLIED CFL NUMBER)			*/
	CFL_max							= 100.0;				/* MAXIMUM CFL NUMBER																*/
	max_dt							= 1.0e-3;				/* LIMITATION OF TIME STEP															*/

	pressure_switch_factor			= 6.0;					/* SET PRESSURE SWITCH FACTOR														*/
	boundary_layer_distance			= 5.0e-6;				/* BOUNDARY LAYER DISTANCE															*/

	relaxation_factor_inviscid		= 1.0;					/* INVISCID JACOBIAN RELAXATION FACTOR FOR STABILITY								*/
	relaxation_factor_viscous		= 1.0;					/* VISCOUS JACOBIAN RELACATION FACTOR FOR STABILITY									*/
	CFL								= 0.0001;				/* CFL NUMBER FOR CURRENT ITERATION													*/
}

PROBLEM_CFD::~PROBLEM_CFD()
{
}

void PROBLEM_CFD::read_CFD(string filename)
{
	int error_code;
	string line;
	string error_message;
	ifstream prob_file;

	prob_file.open(filename.c_str());


	// READ PROBLEM SETUP DATA FILE
	if (prob_file.is_open())
	{
		while (! prob_file.eof())
		{
			getline(prob_file, line);
			data_read::read_line(line);

			data_read::get_data_from_string<int>(line,		"NUMBER OF FLUID:", 0, multi_fluid);

			// NUMERICAL METHOD SETTING - CFD
			data_read::get_data_from_string<int>(line,		"SPARTIAL INTEGRATION METHOD:", 0, SPATIAL_INTEGRATION_METHOD);
			data_read::get_data_from_string<int>(line,		"TIME INTEGRATION METHOD:", 0, TIME_INTEGRATION_METHOD);
			data_read::get_data_from_string<int>(line,		"NUMERICAL ORDER:", 1, NUMERICAL_ORDER);
			data_read::get_data_from_string<int>(line,		"LIMITER METHOD:", 0, LIMITER);

			data_read::get_data_from_string<double>(line,		"RELAXATION FACTOR FOR INVISCID JACOBIAN:", 1.0, relaxation_factor_inviscid);
			data_read::get_data_from_string<double>(line,		"RELAXATION FACTOR FOR VISCOUS JACOBIAN:", 1.0, relaxation_factor_viscous);
		}


		// Read details for SW-FVS method
		if (SPATIAL_INTEGRATION_METHOD == SW_FVS_LAURA || SPATIAL_INTEGRATION_METHOD == SW_FVS_DPLR)
		{
			prob_file.clear();
			prob_file.seekg(0);
			while (! prob_file.eof())
			{
				getline(prob_file, line);
				data_read::read_line(line);

				data_read::get_data_from_string<double>(line,	"PRESSURE SWITCH FACTOR:", 	   6.0, pressure_switch_factor);
				data_read::get_data_from_string<double>(line,	"BOUNDARY LAYER DISTANCE:", 5.0e-6, boundary_layer_distance);
			}
		}


		// Details for implicit methods
		if (TIME_INTEGRATION_METHOD != EXPLICIT)
		{
			prob_file.clear();
			prob_file.seekg(0);
			while (! prob_file.eof())
			{
				getline(prob_file, line);
				data_read::read_line(line);

				data_read::get_data_from_string<int>(line,		"CFL NUMBER INCREASE METHOD:", 0, CFL_INCREASE_METHOD);
				data_read::get_data_from_string<int>(line,		"NUMBER OF ITERATION BEFORE CFL < 1:", 0, iteration_before_1);
				data_read::get_data_from_string<double>(line,	"INITIAL CFL NUMBER:", 0.001, CFL_start);
				data_read::get_data_from_string<double>(line,	"MAXIMUM CFL NUMBER:", 100.0, CFL_max);
				data_read::get_data_from_string<double>(line,	"LIMITATION OF TIME STEP, [SEC]:", 1.0, max_dt);
			}
		}
		else
		{
			prob_file.clear();
			prob_file.seekg(0);
			while (! prob_file.eof())
			{
				getline(prob_file, line);
				data_read::read_line(line);

				data_read::get_data_from_string<double>(line,	"INITIAL CFL NUMBER:", 0.001, CFL_start);
				data_read::get_data_from_string<double>(line,	"LIMITATION OF TIME STEP, [SEC]:", 1.0, max_dt);
			}
		}
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "OP2A_problem_CFD";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A PROBLEM SETUP DATA FILE. PLEASE CHECK THE FILE";
		error.print_message();
	}

	CFL	= CFL_start;

	prob_file.close();
}
