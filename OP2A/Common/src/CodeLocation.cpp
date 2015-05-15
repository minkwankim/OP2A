/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 11, 2015
 *      			Author: Minkwan Kim
 *
 * CodeLocation.hpp
 * 			-  
 *  
 */




#include "Common/include/CodeLocation.hpp"



namespace OP2A {
namespace Common {

/*
 * This class stores the information about a location in the source code
 */
Code_location::Code_location(const char * file, int line, const char * function)
 : m_file(file), m_function(function), m_line (line)
{
}



std::string Code_location::str () const
{
	std::string place (m_file);
	place += ":";
	place += m_line;

	if (!m_function.empty()) // skip if compiler doees not set function
	{
		place += ":";
		place += m_function;
	}

	return place;
}

}
}
