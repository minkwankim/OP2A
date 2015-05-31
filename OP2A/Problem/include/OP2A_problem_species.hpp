/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_problem_species.hpp
 * 			-  
 *  
 */

#ifndef OP2A_PROBLEM_SPECIES_HPP_
#define OP2A_PROBLEM_SPECIES_HPP_


#include <string>
#include <vector>

using namespace std;



/*
 * ======================================
 * 		PROBLEM SETUP for Species DATA
 * ======================================
 */
class PROBLEM_SPECIES_BASIC
{
public:
	int		NS;
	string	species_file;												/* Species data file name			*/

	PROBLEM_SPECIES_BASIC();
	~PROBLEM_SPECIES_BASIC();

	void read_SPECIES_BASIC(string filename);
};




#endif /* OP2A_PROBLEM_SPECIES_HPP_ */
