/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: May 26, 2014
 *      			Author: Minkwan Kim
 *
 * version.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_version.hpp"

////////////////////////////////
/* Member function definitions */
/////////////////////////////////
// Constructor and desctructor
Ver_type::Ver_type(void)
{
	primary		= 0;       /* Primary version number   (i.e. 2 for V2.0) */
	secondary	= 0;     /* Secondary version number (i.e. 0 for V2.0) */
	year		= 0;          /* Most recent date of update */
	month		= 0;
	date		= 0;
}
Ver_type::~Ver_type(void)
{
	primary		= NULL;       /* Primary version number   (i.e. 2 for V2.0) */
	secondary	= NULL;     /* Secondary version number (i.e. 0 for V2.0) */
	year		= NULL;          /* Most recent date of update */
	month		= NULL;
	date		= NULL;
}

// ::M-01
// Assign values
void Ver_type::set_info(int ver_main, int ver_sub, int date_year, int date_month, int date_day, string type_info)
{
	primary		= ver_main;
	secondary	= ver_sub;
	year		= date_year;
	month		= date_month;
	date		= date_day;

	type		= type_info;
}

// ::M-02
// Show version information
void Ver_type::out_info(void)
{
	cout << "OPen source multi-Physics Analyzer (OPPA)    Version "<< primary << "." <<secondary << " [Type]: " << type << endl << endl;
	cout << "Copyright (c) 2013-" << year << "  Min Kwan Kim" << endl <<endl;
	cout << "[IMPORTANT] It is open-source code"<< endl << endl << endl;
}

void Ver_type::out_info(int ver_main, int ver_sub, int date_year, int date_month, int date_day, string type_info)
{
	set_info(ver_main, ver_sub, date_year, date_month, date_day, type_info);
	out_info();
}
