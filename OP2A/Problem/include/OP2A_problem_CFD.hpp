/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_problem_CFD.hpp
 * 			-  
 *  
 */

#ifndef OP2A_PROBLEM_CFD_HPP_
#define OP2A_PROBLEM_CFD_HPP_

#include <string>
#include <vector>


using namespace std;


/*
 * ============================
 * 		PROBLEM SETUP for CFD
 * ============================
 */
class PROBLEM_CFD
{
public:
	int		multi_fluid;
	int		SPATIAL_INTEGRATION_METHOD;				/* SPARTIAL INTEGRATION METHOD		*/
	int		TIME_INTEGRATION_METHOD;				/* TIME INTEGRATION METHOD			*/
	int		NUMERICAL_ORDER;						/* NUMERICAL ORDER					*/
	int		LIMITER;								/* LIMITER METHODS					*/
	int		USE_least_square_method;

	int		CFL_INCREASE_METHOD;					/* CFL NUMBER INCREASE METHOD (IT IS ONLY USED WHEN IMPLICIT METHODS ARE EMPLOYED)	*/
	int		iteration_before_1;						/* NUMBER OF INTERATION BEFORE CFL NUMBER LARGER THAN 1 IS APPLIED					*/
	double	CFL_start;								/* INITIAL CFL NUMBER (FOR EXPLICIT METHOD, IT IS AN APPLIED CFL NUMBER)			*/
	double	CFL_max;								/* MAXIMUM CFL NUMBER																*/
	double	max_dt;									/* LIMITATION OF TIME STEP															*/

	double	pressure_switch_factor;					/* SET PRESSURE SWITCH FACTOR														*/
	double	boundary_layer_distance;				/* BOUNDARY LAYER DISTANCE															*/

	double 	relaxation_factor_inviscid;				/* INVISCID JACOBIAN RELAXATION FACTOR FOR STABILITY								*/
	double 	relaxation_factor_viscous;				/* VISCOUS JACOBIAN RELACATION FACTOR FOR STABILITY									*/
	double 	CFL;									/* CFL NUMBER FOR CURRENT ITERATION													*/

	PROBLEM_CFD();
	~PROBLEM_CFD();

	void read_CFD(string filename);
};



#endif /* OP2A_PROBLEM_CFD_HPP_ */
