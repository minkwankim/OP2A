/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 11, 2014
 *      			Author: Minkwan Kim
 *
 * error_code.hpp
 * 			-  
 *  
 */

#ifndef ERROR_CODE_HPP_
#define ERROR_CODE_HPP_



#include <string>

#define ERROR_CFD_FVS				1

#define	ERROR_MODULE_CFD			0
#define	ERROR_MODULE_DSMC			1
#define	ERROR_MODULE_HYPERSONIC		2
#define	ERROR_MODULE_PLASMA			3

#define ERROR_CODE_NO_FIND			1



using namespace std;

class Error_message_type
{
public:
	int		CPU_number;

	int		module_code;
	string 	module_name;

	int 	primary_code;
	int		secondary_code;

	int		location_primary;
	string	location_primary_name;

	int		location_secondary;
	string	location_secondary_name;

	string 	message;

	Error_message_type();
	~Error_message_type();

	/*
	 * Internal functions
	 */
	void print_message();
	void asg_module_name();
};

void program_error_type1(string error_message);

bool error_check_double(double Value);
bool error_check_int(int Value);
bool error_check_double_pos(double Value);
bool error_check_double_range_bt(double Value, double A);
bool error_check_double_range_st(double Value, double A);







#endif /* ERROR_CODE_HPP_ */
