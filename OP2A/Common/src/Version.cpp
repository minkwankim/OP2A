/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 11, 2015
 *      			Author: Minkwan Kim
 *
 * Version.cpp
 * 			-  
 *  
 */





#include "Common/include/Version.hpp"


namespace OP2A {
namespace Common {

Version::Version(int ver_main, int ver_sub, int date_year, int date_month, int date_day, const char * type_info)
:primary(ver_main), secondary(ver_sub), year(date_year), month(date_month), date(date_day), type(type_info)
{

}

void Version::info(void)
{
	std::cout << "Open-source multi-Physics Phenomena Analyzer (OP2A)	Version "<< primary << "." <<secondary << " [Type]: " << type << std::endl << std::endl;
	std::cout << "Copyright (c) 2013-" << year << "  Min Kwan Kim" << std::endl << std::endl;
	std::cout << "[IMPORTANT] It is open-source code"<< std::endl << std::endl << std::endl;
}





}
}
