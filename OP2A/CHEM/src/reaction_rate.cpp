/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Feb 11, 2015
 *      			Author: Minkwan Kim
 *
 * reaction_rate.cpp
 * 			-  
 *  
 */



#include "../include/reaction.hpp"
#include "../include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"


REACTION_RATE::REACTION_RATE()
{
	method	= 0;
	rate_coeff.resize(4);
	temperature.resize(4);
}

REACTION_RATE::~REACTION_RATE()
{

}


double REACTION_RATE::calculate_rate_coefficient(double T)
{
	double kf = 0.0;

	kf	= rate_coeff[0]	* pow(T, rate_coeff[1])	* exp(-rate_coeff[2] * pow(T, rate_coeff[3]));
	return (kf);
}

double REACTION_RATE::calculate_temperature(vector<double> &T)
{
	double sum_coeff = 0.0;
	for (int i = 0; i <= 3; i++)
	{
		if (T[i] > 0.0) sum_coeff += temperature[i];
	}
	for (int i = 0; i <= 3; i++)	temperature[i] /= sum_coeff;

	double Tp	= pow(T[0], temperature[0]);
	for (int i = 1; i <= 3; i++)
	{
		if (T[i] > 0.0)	Tp = Tp * pow(T[i], temperature[i]);
	}

	return (Tp);
}
