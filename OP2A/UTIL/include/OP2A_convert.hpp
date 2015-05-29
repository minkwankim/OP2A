/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_convert.hpp
 * 			-  
 *  
 */

#ifndef OP2A_CONVERT_HPP_
#define OP2A_CONVERT_HPP_



#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <cmath>
#include <memory>
#include <limits>


#ifdef MPIF
	#include <mpi.h>
#endif

using namespace std;



string convertInt(int number);
string convertDouble(double number);
string convert_YES_NO(int i);



#endif /* OP2A_CONVERT_HPP_ */
