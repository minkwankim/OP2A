/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_problem_viscous.hpp
 * 			-  
 *  
 */

#ifndef OP2A_PROBLEM_VISCOUS_HPP_
#define OP2A_PROBLEM_VISCOUS_HPP_



#include <string>
#include <vector>

using namespace std;


/*
 * ======================================
 * 		PROBLEM SETUP for VISCOUS FLOW
 * ======================================
 */
class PROBLEM_VISCOUS
{
public:
	int	model_viscosity;
	int	model_conductivity;
	int model_mixing_rule;

	int 	viscous_relaxation_na;
	int 	viscous_relaxation_nb;
	double 	viscous_relaxation_min;
	double 	viscous_relaxation_max;

	double 	Le;

	bool adiabatic;
	bool catalytic;
	bool radiative;
	bool use_emissivity;

	PROBLEM_VISCOUS();
	~PROBLEM_VISCOUS();

	void read_VISCOUS(string filename);
};


#endif /* OP2A_PROBLEM_VISCOUS_HPP_ */
