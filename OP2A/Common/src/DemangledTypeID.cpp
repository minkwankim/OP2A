/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 26, 2015
 *      			Author: Minkwan Kim
 *
 * DemangledTypeID.cpp
 * 			-  
 *  
 */

#include "Common/include/OP2A.hpp"
#include "Common/include/DemangledTypeID.hpp"


namespace OP2A{
namespace Common{

std::string demangle(const char* type)
{
	int status = 0;
	char* r = 0;

	std::string ret_value;
	if ( (r == 0) || (status != 0) )	ret_value = std::string(type);
	else 								ret_value = std::string(r);

	free(r);

	return ret_value;
}



}
}

