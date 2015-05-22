/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Feb 11, 2015
 *      			Author: Minkwan Kim
 *
 * reaction_keq.cpp
 * 			-  
 *  
 */



#include "../include/reaction.hpp"
#include "../include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../../MATRIX/include/math_misc.hpp"


REACTION_KEQ::REACTION_KEQ()
{
	num_data = 0;
	model	 = 0;
}

REACTION_KEQ::~REACTION_KEQ()
{

}



double REACTION_KEQ::calculate_logKeq(double T, double n_mix)
{
	double a1	= interpolate2 (n, A1, num_data, n_mix);
	double a2	= interpolate2 (n, A2, num_data, n_mix);
	double a3	= interpolate2 (n, A3, num_data, n_mix);
	double a4	= interpolate2 (n, A4, num_data, n_mix);
	double a5	= interpolate2 (n, A5, num_data, n_mix);

	double log_keq;
	double Z	= 10000.0/T;
	switch (model)
	{
	case PARK85:
		log_keq	= a1 + a2*Z + a3*Z*Z + a4*pow(Z, 3.0) + a5*pow(Z, 4.0);
		break;

	case PARK90:
		log_keq	= a1/Z + a2 + a3*log(Z) + a4*Z + a5*Z*Z;
		break;

	case PARK94:
		log_keq	= a1 + a2*log(Z) + a3*Z + a4*pow(Z, 2.0) + a5*pow(Z, 3.0);
		break;
	}

	if (log_keq != log_keq || fabs(log_keq) == numeric_limits<double>::infinity())
	{
		Error_message_type	error_message;

		error_message.module_name	="CHEMISTRY-REACTION";
		error_message.location_primary_name	= "N/A";
		error_message.location_secondary_name	= "N/A";
		error_message.message	= " Problem in the calculation of equilibrium constant!";
		error_message.print_message();
	}

	return (log_keq);
}







double REACTION_KEQ::calculate_dKeq_dT_over_Keq(double T, double n_mix)
{
	double a1	= interpolate2 (n, A1, num_data, n_mix);
	double a2	= interpolate2 (n, A2, num_data, n_mix);
	double a3	= interpolate2 (n, A3, num_data, n_mix);
	double a4	= interpolate2 (n, A4, num_data, n_mix);
	double a5	= interpolate2 (n, A5, num_data, n_mix);

	double dKeq_dT_over_Keq;
	double Z	= 10000.0/T;
	switch (model)
	{
	case PARK85:
		dKeq_dT_over_Keq	= -(a2*Z + 2.0*a3*Z*Z + 3.0*a4*pow(Z, 3.0) + 4.0*a5*pow(Z, 4.0)) / T;
		break;

	case PARK90:
		dKeq_dT_over_Keq	= -(a1/Z + a3 + a4*Z + 2.0*a5*Z*Z) / T;
		break;

	case PARK94:
		dKeq_dT_over_Keq	= -(a2 + a3*Z + 2.0*a4*Z*Z + 3.0*a5*pow(Z, 3.0)) / T;
		break;
	}

	if (dKeq_dT_over_Keq != dKeq_dT_over_Keq || fabs(dKeq_dT_over_Keq) == numeric_limits<double>::infinity())
	{
		Error_message_type	error_message;

		error_message.module_name	="CHEMISTRY-REACTION";
		error_message.location_primary_name	= "N/A";
		error_message.location_secondary_name	= "N/A";
		error_message.message	= " Problem in the calculation of dKeq_dT_over_Keq!";
		error_message.print_message();
	}

	return (dKeq_dT_over_Keq);
}



