/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_WIG_Problem.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_WIG_Problem.hpp"


PROBLEM_SPECIES::PROBLEM_SPECIES()
{
	NS_n	= 0;
	NS_i	= 0;
}

PROBLEM_SPECIES::~PROBLEM_SPECIES()
{

}


void PROBLEM_SPECIES::read_SPECIES(string file_name)
{
	read_SPECIES_BASIC(file_name);
}






PROBLEM_CLASS_PLASMA_ver1::PROBLEM_CLASS_PLASMA_ver1()
{

}


PROBLEM_CLASS_PLASMA_ver1::~PROBLEM_CLASS_PLASMA_ver1()
{

}



void PROBLEM_CLASS_PLASMA_ver1::read(string file_name)
{
	read_COMMON(file_name);
	read_CFD(file_name);
	read_SPECIES(file_name);
	read_NONEQ(file_name);
	read_VISCOUS(file_name);
	read_WALLCOND(file_name);

	IC.read_IC(file_name);


	// Processing Initial conditions
	for (int n	= 0; n <= IC.NIC-1; n++)
	{
		IC.rho[n]	= 0.0;
		for (int s = 0; s <= NS-1; s++)
		{
			IC.rho[n]	+= IC.rho_s[n][s];
		}
	}

}



// M::02 - Adjust problem setting
void PROBLEM_CLASS_PLASMA_ver1::adjust (SPECIES_DATA_BASIC species)
{
	vector	<vector <double> >		rho_s_temp = IC.rho_s;
	vector	<vector <double> >		Ys_temp   = Ys;

	for (int i = 0; i <= IC.NIC-1; i++)
	{
		for (int s = 0; s <= MAX_NS_PROBLEM-1; s++)	IC.rho_s[i][s] = 0.0;

		IC.rho[i] = 0.0;
		for (int s = 0; s <= species.NS-1; s++)
		{
			int s_ID = species.whereis[s];
			IC.rho_s[i][s]	= rho_s_temp[i][s_ID];
			IC.rho[i] += IC.rho_s[i][s];
		}
	}


	for (int i = 0; i <= NWCOND-1; i++)
	{
		for (int s = 0; s <= MAX_NS_PROBLEM-1; s++)	Ys[i][s] = 0.0;

		double Ys_tot = 0.0;
		for (int s = 0; s <= species.NS-1; s++)
		{
			int s_ID 	= species.whereis[s];
			Ys[i][s]	= Ys_temp[i][s_ID];
			Ys_tot 		+= Ys[i][s];
		}

		for (int s = 0; s <= species.NS-1; s++)	Ys[i][s] = Ys[i][s] / Ys_tot;
	}
}
