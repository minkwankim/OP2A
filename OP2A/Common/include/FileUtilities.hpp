/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 26, 2015
 *      			Author: Minkwan Kim
 *
 * FileUtilities.hpp
 * 			-  
 *  
 */
#ifndef FILEUTILITIES_HPP_
#define FILEUTILITIES_HPP_

#include <fstream>

namespace OP2A{

namespace Common{

bool fileExist(const char *filename);
bool fileExist(const std::string &filename);



}
}



#endif /* FILEUTILITIES_HPP_ */
