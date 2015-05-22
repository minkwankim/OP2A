/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 7, 2015
 *      			Author: Minkwan Kim
 *
 * reaction_rate_coeff.cpp
 * 			-  
 *  
 */


#include <omp.h>

#include "../include/reaction.hpp"
#include "../include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"




/*
 * Constructor
 */
 REACTION_RATE_COEFF::REACTION_RATE_COEFF()
 {
	 method	= 0;				// CALCULATION METHOD

	 for (int i = 0; i <= 3; i++)
	 {
		 kf_coeff[i]			= 0.0;;
		 temperature_coeff_f[i]	= 0.0;;
	 }

	 for (int i = 0; i <= 3; i++)
	 {
		 kb_coeff[i]			= 0.0;
		 temperature_coeff_b[i]	= 0.0;;
	 }
 }

 REACTION_RATE_COEFF::~REACTION_RATE_COEFF()
 {

 }




 /*
  * Calculate Forward reaction rare coefficient
  */
 double cal_kf(REACTION_RATE_COEFF reaction_k, double Tf)
 {
	 double kf = 0.0;
	 kf	= reaction_k.kf_coeff[0]	* pow(Tf, reaction_k.kf_coeff[1])	* exp(-reaction_k.kf_coeff[2] * pow(Tf, reaction_k.kf_coeff[3]));

	if(kf != kf	||  kf == numeric_limits<double>::infinity())
	{
		Error_message_type	error_message;
		error_message.module_name	="Chemistry-Reaction_rate_coeff";
		error_message.location_primary_name	= "N/A";
		error_message.location_secondary_name	= "N/A";
		error_message.message	= " Problem in forward reaction rate!";
		error_message.print_message();
	}

	return (kf);
 }


 /*
  * Calculate Backward reaction rare coefficient
  */
double cal_kb(REACTION_RATE_COEFF reaction_k, double n_mix, double Tb)
{
	double kb = 0.0;
	switch(reaction_k.method)
	{
	case 0:
		kb	= reaction_k.kb_coeff[0]	* pow(Tb, reaction_k.kb_coeff[1])	* exp(-reaction_k.kb_coeff[2] * pow(Tb, reaction_k.kb_coeff[3]));
		break;
	case 1:
		double log_Keq	= reaction_k.Keq_data.calculate_logKeq(Tb, n_mix);
		double kfb;
		kfb	= reaction_k.kf_coeff[0]	* pow(Tb, reaction_k.kf_coeff[1])	* exp(-reaction_k.kf_coeff[2] * pow(Tb, reaction_k.kf_coeff[3]));

		if (kfb != 0.0)
		{
			double log_kfb	= log(kfb);
			double log_kb	= log_kfb - log_Keq;

			if (log_kb != log_kb  || fabs(log_kb) == numeric_limits<double>::infinity())
			{
				Error_message_type	error_message;
				error_message.module_name	="Chemistry-Reaction_rate_coeff";
				error_message.location_primary_name	= "N/A";
				error_message.location_secondary_name	= "N/A";
				error_message.message	= " Problem in backward reaction rate!";
				error_message.print_message();
			}
			kb	= exp(log_kb);
		}
		else
		{
			kb = 0.0;
		}
	}

	kb = 0.0;

	return(kb);
}




/*
 * Calculate reaction temperatures
 */
void calculate_reaction_temperature(REACTION_RATE_COEFF reaction_k, vector<double> &T, double &Tf, double &Tb)
{
	kmp_set_defaults("KMP_AFFINITY = scatter");

	 vector<double> temp_f(4,0);
	 vector<double> temp_b(4,0);

#pragma omp parallel for
	 for (int i = 0; i <= 3; i++)	temp_f[i]	= pow(T[i], temp_f[i]);

#pragma omp parallel for
	 for (int i = 0; i <= 3; i++)	temp_b[i]	= pow(T[i], temp_b[i]);

	 Tf	= temp_f[0] * temp_f[1] * temp_f[2] * temp_f[3];
	 Tb	= temp_b[0] * temp_b[1] * temp_b[2] * temp_b[3];
}





