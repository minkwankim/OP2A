/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 26, 2015
 *      			Author: Minkwan Kim
 *
 * FileUtilities.cpp
 * 			-  
 *  
 */



#include "Common/include/FileUtilities.hpp"

namespace OP2A{
namespace Common{

bool fileExist(const char *filename)
{
	std::ifstream ifile(filename);
	return ifile;
}

bool fileExist(const std::string& filename)
{
  std::ifstream ifile(filename.c_str());
  return ifile;
}

}
}
