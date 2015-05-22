/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_problem_NONEW.hpp
 * 			-  
 *  
 */

#ifndef OP2A_PROBLEM_NONEQ_HPP_
#define OP2A_PROBLEM_NONEQ_HPP_


#include <string>
#include <vector>

#define	MAX_NS_PROBLEM	20

using namespace std;


/*
 * ============================
 * 		PROBLEM SETUP for NONEQ
 * ============================
 */
class PROBLEM_NONEQ
{
public:
	int	NER;
	int	NEV;
	int	NEE;

	int	ID_Tr;
	int	ID_Tv;
	int	ID_Te;

	bool 	NONEQ_Chem;
	string	reaction_file;

	PROBLEM_NONEQ();
	~PROBLEM_NONEQ();

	void read_NONEQ(string filename);
};





#endif /* OP2A_PROBLEM_NONEW_HPP_ */
