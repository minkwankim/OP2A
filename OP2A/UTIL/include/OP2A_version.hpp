/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: May 26, 2014
 *      			Author: Minkwan Kim
 *
 * version.hpp
 * 			-  
 *  
 */

#ifndef OPPA_VERSION_HPP_
#define OPPA_VERSION_HPP_


#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <cmath>
#include <memory>




#ifdef MPI
	#include <mpi.h>
#endif

using namespace std;


/* Version information */
class Ver_type
{
public:
	int primary;       /* Primary version number   (i.e. 2 for V2.0) */
	int secondary;     /* Secondary version number (i.e. 0 for V2.0) */
	int year;          /* Most recent date of update */
	int month;
	int date;

	string type;     /* "Type" of version */

	// Constructor and desctructor
	Ver_type(void);
	~Ver_type(void);

	// ::M-01
	// Assign values
	void set_info(int ver_main, int ver_sub, int date_year, int date_month, int date_day, string type_info);

	// ::M-02
	void out_info(void);
	void out_info(int ver_main, int ver_sub, int date_year, int date_month, int date_day, string type_info);
};

#endif /* VERSION_HPP_ */
