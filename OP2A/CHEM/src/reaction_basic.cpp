/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: Nov 21, 2014
 *      			Author: Minkwan Kim
 *
 * reaction.cpp
 * 			-  
 *  
 */

#include "../include/reaction.hpp"
#include "../include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"



REACTION_BASIC::REACTION_BASIC()
{
	ID			= -1;
	type		= -1;				// Reaction type
	Reactant_num	= 0;		// Number of reactant
	Product_num		= 0;
}

REACTION_BASIC::~REACTION_BASIC()
{

}





REACTION_BASIC_ver2::REACTION_BASIC_ver2()
{
	ID			= -1;
	type		= -1;				// Reaction type
}

REACTION_BASIC_ver2::~REACTION_BASIC_ver2()
{

}








/*

Reaction_plasma::Reaction_plasma()
{
	type 			= -1;
	Reactant_num	= 0;
	Product_num		= 0;

	// Reaction coefficient
	method_kb	= 0;
	Af	= 0.0;
	Bf	= 0.0;
	Cf	= 0.0;
	Df	= 0.0;

	Ab	= 0.0;
	Bb	= 0.0;
	Cb	= 0.0;
	Db	= 0.0;

	af	= 0.0;
	bf	= 0.0;
	cf	= 0.0;
	df	= 0.0;
	ef	= 0.0;

	ab	= 0.0;
	bb	= 0.0;
	cb	= 0.0;
	db	= 0.0;
	eb	= 0.0;

	aK	= 0.0;
	bK	= 0.0;
	cK	= 0.0;
	dK	= 0.0;
	eK	= 0.0;
}


Reaction_plasma::~Reaction_plasma()
{

}





double Reaction_plasma::kf(double T)
{
	double reaction_rate = 0.0;
	reaction_rate	= Af * pow(T, Bf) * exp(-Cf*pow(T, Df));
	return(reaction_rate);
}

double Reaction_plasma::kb(double T)
{
	double reaction_rate = 0.0;
	reaction_rate	= Ab * pow(T, Bb) * exp(-Cb*pow(T, Db));
	return(reaction_rate);
}


double Reaction_plasma::Tp_forward(vector <vector <double> > &T, SPECIES_DATA_BASIC &species_entire)
{
	double Tm;
	vector<double> Tmean(5, 1.0);


	int ns	= 0;
	for (int s = 0; s <= Reactant_num-1; s++)
	{
		int s_ID	= Reactant_id[s];
		int s_ptr	= species_entire.whereis[s_ID];

		if ((*species_entire.data_entire)[s_ptr].basic_data.type != ELECTRON)
		{
			for (int k = 0; k <= 3; k++)	Tmean[k] *= T[s_ptr][k];
			ns++;
		}
		else
		{
			Tmean[4]	= T[s_ptr][0];
		}
	}

	for (int k = 0; k <= 3; k++)	Tmean[k] = pow(Tmean[k], 1.0/ns);

	Tm	= pow(Tmean[0], af) * pow(Tmean[1], bf) * pow(Tmean[2], cf) * pow(Tmean[3], df) * pow(Tmean[4], ef);
	return (Tm);
}


double Reaction_plasma::Tp_backward(vector <vector <double> > &T, SPECIES_DATA_BASIC &species_entire)
{
	double Tm;
	vector<double> Tmean(5, 1.0);


	int ns	= 0;
	for (int s = 0; s <= Product_num-1; s++)
	{
		int s_ID	= Product_id[s];
		int s_ptr	= species_entire.whereis[s_ID];

		if ((*species_entire.data_entire)[s_ptr].basic_data.type != ELECTRON)
		{
			for (int k = 0; k <= 3; k++)	Tmean[k] *= T[s_ptr][k];
			ns++;
		}
		else
		{
			Tmean[4]	= T[s_ptr][0];
		}
	}

	for (int k = 0; k <= 3; k++)	Tmean[k] = pow(Tmean[k], 1.0/ns);

	Tm	= pow(Tmean[0], ab) * pow(Tmean[1], bb) * pow(Tmean[2], cb) * pow(Tmean[3], db) * pow(Tmean[4], eb);
	return (Tm);
}

double Reaction_plasma::Tp_Keq(vector <vector <double> > &T, SPECIES_DATA_BASIC &species_entire)
{
	double Tm;
	vector<double> Tmean(5, 1.0);


	int ns	= 0;
	for (int s = 0; s <= Reactant_num-1; s++)
	{
		int s_ID	= Reactant_id[s];
		int s_ptr	= species_entire.whereis[s_ID];

		if ((*species_entire.data_entire)[s_ptr].basic_data.type != ELECTRON)
		{
			for (int k = 0; k <= 3; k++)	Tmean[k] *= T[s_ptr][k];
			ns++;
		}
		else
		{
			Tmean[4]	= T[s_ptr][0];
		}
	}

	for (int s = 0; s <= Product_num-1; s++)
	{
		int s_ID	= Product_id[s];
		int s_ptr	= species_entire.whereis[s_ID];

		if ((*species_entire.data_entire)[s_ptr].basic_data.type != ELECTRON)
		{
			for (int k = 0; k <= 3; k++)	Tmean[k] *= T[s_ptr][k];
			ns++;
		}
		else
		{
			Tmean[4]	= T[s_ptr][0];
		}
	}


	for (int k = 0; k <= 3; k++)	Tmean[k] = pow(Tmean[k], 1.0/ns);

	Tm	= pow(Tmean[0], aK) * pow(Tmean[1], bK) * pow(Tmean[2], cK) * pow(Tmean[3], dK) * pow(Tmean[4], eK);
	return (Tm);
}
*/
