/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 20, 2015
 *      			Author: Minkwan Kim
 *
 * LogFile.cpp
 * 			-  
 *  
 */


#include "Common/include/LogFile.hpp"

namespace OP2A{


// Constructo9r
Logger::Logger()
{


}

Logger::~Logger()
{

}



Logger& Logger::getInstance()
{
	static Logger logger;
	return logger;
}



}
