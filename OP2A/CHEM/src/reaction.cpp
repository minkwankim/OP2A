/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Feb 12, 2015
 *      			Author: Minkwan Kim
 *
 * reaction.cpp
 * 			-  
 *  
 */




#include "../include/reaction.hpp"
#include "../include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"


REACTION::REACTION()
{
	Tref	= 0.0;
}

REACTION::~REACTION()
{

}




double REACTION::cal_kf(vector< vector<double> > &T_All, double n_mix)
{
	vector<double>	T(4, 1.0);

	for (int k = 0; k <= 3; k++)
	{
		for (int i = 0; i <= Reactant_num-1; i++)
		{
			int s	= Reactant_id[i];

			if (T_All[s][k] > 0.0)	T[k]	*= T_All[s][k];
		}

		T[k]	= pow(T[k], 1.0/Reactant_num);
	}

	double T_kf;
	T_kf	= kf.calculate_temperature(T);


	double k_forward	= kf.calculate_rate_coefficient(T_kf);
	return (k_forward);
}



double REACTION::cal_kb(vector< vector<double> > &T_All, double n_mix)
{
	vector<double>	T(4, 1.0);

	for (int k = 0; k <= 3; k++)
	{
		for (int i = 0; i <= Product_num-1; i++)
		{
			int s	= Product_id[i];

			if (T_All[s][k] > 0.0)	T[k]	*= T_All[s][k];
		}

		T[k]	= pow(T[k], 1.0/Product_num);
	}


	double T_kb;
	T_kb	= kb.calculate_temperature(T);

	double k_backward = 0.0;

	if (kb.method == 0)
	{
		k_backward = kb.calculate_rate_coefficient(T_kb);
	}
	else
	{
		double log_keq 	= Keq.calculate_logKeq(T_kb, n_mix);
		double log_kf	= kf.calculate_rate_coefficient(T_kb);
		log_kf	= log(log_kf);

		double log_kb	= log_kf - log_keq;

		if (log_kb != log_kb  || fabs(log_keq) == numeric_limits<double>::infinity())
		{
			Error_message_type	error_message;

			error_message.module_name	="CHEMISTRY-REACTION - BACKWARD kb";
			error_message.location_primary_name	= "N/A";
			error_message.location_secondary_name	= "N/A";
			error_message.message	= " Problem in the calculation of backward reaction rate!";
			error_message.print_message();
		}

		k_backward	= exp(log_kb);
	}

	return (k_backward);
}









REACTION_SINGLE::REACTION_SINGLE()
{
	Tref	= 0.0;
}

REACTION_SINGLE::~REACTION_SINGLE()
{

}
