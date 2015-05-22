/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_problem_IC.hpp
 * 			-  
 *  
 */

#ifndef OP2A_PROBLEM_IC_HPP_
#define OP2A_PROBLEM_IC_HPP_



#include <string>
#include <vector>

using namespace std;


/*
 * ======================================
 * PROBLEM SETUP for INITIAL CONDITONS
 * ======================================
 */
class PROBLEM_IC_BASIC
{
public:

	// Initial Conditions and Method
	int					INITIALIZE_METHOD;						/* METHOD FOR INITIALIZATION => 0: USING INFLOW CONDITIONS / 1: USING AMBIENT CONDITIONS */
	unsigned int		NIC;									/* Number of conditions */
	vector <int>		whereis;

	vector	<vector <double> >		rho_s;
	vector	<double>				rho;
	vector	<vector <double> >		v;
	vector	<double>		T;
	vector	<double>		Tr;
	vector	<double> 		Tv;
	vector	<double>		Te;

	PROBLEM_IC_BASIC();
	~PROBLEM_IC_BASIC();

	// Internal functions
	void read_IC(string filename);
};



class PROBLEM_IC_CFD : public PROBLEM_IC_BASIC
{
public:
	vector <double> E;
	vector <double> Er;
	vector <double> Ev;
	vector <double> Ee;

	PROBLEM_IC_CFD();
	~PROBLEM_IC_CFD();

	void assign_size();
};




#endif /* OP2A_PROBLEM_IC_HPP_ */
