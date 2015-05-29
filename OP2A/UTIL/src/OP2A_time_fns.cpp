/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: Oct 24, 2014
 *      			Author: Minkwan Kim
 *
 * misc_fns.cpp
 * 			-  
 *  
 */

#include <vector>
#include <sys/time.h>

using namespace std;

double dtime()
{
	double tseconds	= 0.0;
	struct timeval	mytime;

	gettimeofday(&mytime, (struct timezone*)0);

	tseconds	= (double)(mytime.tv_sec	+ mytime.tv_usec*1.0e-6);
	return	(tseconds);
}

