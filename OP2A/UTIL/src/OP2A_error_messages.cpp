/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: Dec 19, 2014
 *      			Author: Minkwan Kim
 *
 * error_fns.cpp
 * 			-  
 *  
 */

#include <iostream>
#include <cstdlib>
#include <limits>
#include <math.h>

#include "../include/OP2A_error_message.hpp"

using namespace std;

Error_message_type::Error_message_type()
{
	CPU_number		= -1;
	module_code		= -1;
	primary_code	= -1;
	secondary_code	= -1;

	location_primary	= -1;
	location_secondary	= -1;
}

Error_message_type::~Error_message_type()
{

}




/*
 * Internal functions
 */
void Error_message_type::print_message()
{
	asg_module_name();

	cout << "==========================================================" << endl;
	cout << "[!!!! ERROR AT PROCESSOR NUMBER: " << CPU_number << " !!!!]" << endl;
	cout << "==========================================================" << endl;
	cout << "[Error module]: \t" << module_name << endl;
	cout << "[Primary Error Code]: \t" <<  primary_code << endl;
	cout << "[Secondary Error Code]: \t" <<  secondary_code << endl;
	cout <<  endl;
	cout << "[Error Location]" << endl;
	cout << " \t \t" << location_primary_name << ": \t" << location_primary << endl;
	cout << " \t \t" << location_secondary_name << ": \t" << location_secondary << endl;
	cout << endl;
	cout << "[Error NOTE]:" << endl;
	cout << message << endl;
	cout << endl;
	cout << " ->SIMULATION IS STOPPED" << endl;

#ifdef MPI
	MPI_Abort(MPI_COMM_WORLD, 0);
#endif

	abort();
}


void Error_message_type::asg_module_name()
{
	switch (module_code)
	{
	case ERROR_MODULE_CFD:
		module_name = "CFD";
		break;

	case ERROR_MODULE_DSMC:
		module_name = "DSMC";
		break;

	case ERROR_MODULE_HYPERSONIC:
		module_name = "HYPERSONIC";
		break;

	case ERROR_MODULE_PLASMA:
		module_name = "PLASMA";
		break;
	}
}

// Legacy Error function
void program_error_type1(string error_message)
{
	int myrank;
	myrank = 0;

#ifdef MPIF
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif

	cout << "[!!!! ERROR AT PROCESSOR #[" << myrank << " !!!!]" << endl;
	cout << "[ERROR MESSAGE]: " <<  error_message << endl;
	cout << "SIMULATION IS STOPPED: PLEASE CHECK ERROR MESSAGE TO SOLVE THE PROBLEM!!" << endl;

#ifdef MPIF
	MPI_Abort(MPI_COMM_WORLD, 0);
#endif
	abort();
}

/*
 * ===============================================================================
 * 	Functions related to error check
 * ===============================================================================
 */
bool error_check_double(double Value)
{
	bool flag = false;
	if (Value	!= Value || fabs(Value) == numeric_limits<double>::infinity())	flag = true;

	return (flag);
}

bool error_check_int(int Value)
{
	bool flag = false;
	if (Value	!= Value || abs(Value) == numeric_limits<int>::infinity())	flag = true;

	return (flag);
}

bool error_check_double_pos(double Value)
{
	bool flag = false;
	if (Value	!= Value || Value < 0.0 || fabs(Value) == numeric_limits<double>::infinity())	flag = true;

	return (flag);
}

bool error_check_double_range_bt(double Value, double A)
{
	bool flag = false;
	if (Value	!= Value || Value < A || fabs(Value) == numeric_limits<double>::infinity())	flag = true;

	return (flag);
}

bool error_check_double_range_st(double Value, double A)
{
	bool flag = false;
	if (Value	!= Value || Value > A || fabs(Value) == numeric_limits<double>::infinity())	flag = true;

	return (flag);
}
