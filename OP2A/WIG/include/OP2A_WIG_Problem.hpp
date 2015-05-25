/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OPPA_WIG_Problem.hpp
 * 			-  
 *  
 */

#ifndef OPPA_WIG_PROBLEM_HPP_
#define OPPA_WIG_PROBLEM_HPP_

#include "Problem/include/OP2A_setup_problem.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"




class PROBLEM_SPECIES : public PROBLEM_SPECIES_BASIC
{
public:
	int NS_n;
	int NS_i;

	PROBLEM_SPECIES();
	~PROBLEM_SPECIES();

	void read_SPECIES(string file_name);

};



class PROBLEM_CLASS_PLASMA_ver1 : public PROBLEM_COMMON, public PROBLEM_CFD, public PROBLEM_SPECIES, public PROBLEM_NONEQ, public PROBLEM_VISCOUS, public PROBLEM_WALLCOND
{
public:
	PROBLEM_IC_CFD	IC;

	PROBLEM_CLASS_PLASMA_ver1();
	~PROBLEM_CLASS_PLASMA_ver1();


	// M::01 - Read program input file
	void read(string file_name);

	// M::02 - Adjust problem setting
	void adjust (SPECIES_DATA_BASIC species);
};


#define OP2A_PROBLEM	PROBLEM_CLASS_PLASMA_ver1



#endif /* OPPA_WIG_PROBLEM_HPP_ */
