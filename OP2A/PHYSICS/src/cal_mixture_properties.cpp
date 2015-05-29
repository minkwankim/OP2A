/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * mixture_properties.cpp
 * 			-  
 *  
 */

#include "../include/constant_text.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../include/mixture_properties.hpp"


double cal_mixture_properties_ver1(vector<double>Q, int NS, vector<double> Yi, vector<double> Xi, int mode)
{
	double Q_mix	= 0.0;

	switch (mode)
	{
	case DENSITY:
		for (int s = 0; s <= NS-1; s++)	Q_mix	+= Q[s];
		break;

	case PRESSURE:
		for (int s = 0; s <= NS-1; s++)	Q_mix	+= Q[s];
		break;

	case SPECIFIC_HEAT:
		for (int s = 0; s <= NS-1; s++)	Q_mix	+= Yi[s] * Q[s];
		break;

	case GAS_CONSTANT:
		for (int s = 0; s <= NS-1; s++)	Q_mix	+= (Xi[s] / Q[s]);
		Q_mix	= 1.0 / Q_mix;
		break;

	case MOLECULAR_WEIGHT:
		for (int s = 0; s <= NS-1; s++)	Q_mix	+= Xi[s] * Q[s];
		break;
	}

	return(Q_mix);
}















/*
 * ==============================================================
 * Mixing rules for Transport properties
 * ==============================================================
 */
double mixing_Wilke_main_algorithm(vector<double> Qs, vector<double> Xs, vector<double> phi_s, unsigned int NS)
{
	double Q_mix	= 0.0;

	for (int s = 0; s <= NS-1; s++)	Q_mix	+= (Xs[s] * Qs[s]) / phi_s[s];

	return (Q_mix);
}

void mixing_Wilke_phi_ver1(vector<double> mu_s, vector<double> Xs, vector<double> Ms, vector<double> &phi_s, unsigned int NS)
{
	double temp1;
	double temp2;

	phi_s.resize(NS);

	for (int s = 0; s <= NS-1; s++)
	{
		phi_s[s]	= 0.0;

		for (int r = 0; r <= NS-1; r++)
		{
			temp1	= 1.0 + sqrt(mu_s[s]/mu_s[r]) * pow(Ms[r]/Ms[s], 0.25);
			temp2	= 8.0 * ( 1.0 + Ms[s]/Ms[r]);

			phi_s[s]	+= Xs[r] * (temp1*temp1) / sqrt(temp2);
		}


		// Error Check
		if (error_check_double_pos(phi_s[s] == true))
		{
			Error_message_type	error_message;
			error_message.module_name	="Physics-Mixing Rule: Wilke";
			error_message.location_primary_name	= "Species";
			error_message.location_primary		= s;
			error_message.location_secondary_name	= "NONE";
			error_message.location_secondary		= 0;
			error_message.message	= " Problem in the calculatin of phi_s in Wilke rule";
			error_message.print_message();
		}
	}
}


double mixing_Wilke(vector<double> Qs, vector<double> mu_s, vector<double> Xs, vector<double> Ms, unsigned int NS)
{
	double Q_mix	= 0.0;
	vector <double>	phi_s(NS, 0.0);


	mixing_Wilke_phi(mu_s, Xs, Ms, phi_s, NS);
	Q_mix	= mixing_Wilke_main_algorithm(Qs, Xs, phi_s, NS);

	return (Q_mix);
}
