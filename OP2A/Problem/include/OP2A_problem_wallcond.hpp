/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 13, 2015
 *      			Author: Minkwan Kim
 *
 * problem_class.hpp
 * 			-  
 *  
 */

#ifndef PROBLEM_CLASS_HPP_
#define PROBLEM_CLASS_HPP_

#include <string>
#include <vector>

using namespace std;




/*
 * ======================================
 * 		PROBLEM SETUP for VISCOUS FLOW
 * ======================================
 */
class PROBLEM_WALLCOND
{
public:
	unsigned int				NWCOND;
	vector < vector <double> >	Ys;
	vector < vector <double> >	Tw;

	vector < double> 			emissivity_Tref;
	vector < vector <double> > 	emissivity_below;
	vector < vector <double> > 	emissivity_above;

	PROBLEM_WALLCOND();
	~PROBLEM_WALLCOND();

	void read_WALLCOND(string filename);
};





















/*
 * PROBLEM SETUP for Multi-fluid
 */
class PROBLEM_MULTI_FLUID
{
public:

};


#endif /* PROBLEM_CLASS_HPP_ */
