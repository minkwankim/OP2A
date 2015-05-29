/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_convert.cpp
 * 			-  
 *  
 */


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


string convertInt(int number)
{
   stringstream ss;		//create a stringstream
   ss << number;		//add number to the stream
   return ss.str();		//return a string with the contents of the stream
}

string convertDouble(double number)
{
   stringstream ss;		//create a stringstream
   ss << number;		//add number to the stream
   return ss.str();		//return a string with the contents of the stream
}

string convert_YES_NO(int i)
{
	string ss;

	if (i == 1) ss = "YES";
	else if (i == 0) ss = "NO";

	return ss;
}


