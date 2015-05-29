/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Variable_temp_calculate.cpp
 * 			-  
 *  
 */



#include "../include/OP2A_CFD_Variable_change.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"


/*
 * Variable changes
 */
fn_Evib::fn_Evib()
{
	NS	= 0;
	Eve	= 0;
}

fn_Evib::~fn_Evib()
{

}

void fn_Evib::create(SPECIES_DATA_BASIC &species_data, vector<double> &rho, int ns, double E)
{
	species	= species_data;
	rhos	= rho;
	NS		= ns;
	Eve		= E;
}

double fn_Evib::f(double T)
{
	double fn	= 0.0;

	for (int s = 0; s <= NS-1; s++)
	{
		int ss = species.whereis[s];
		fn += rhos[s]* (*species.data_entire)[ss].e_VEE(T);
	}

	fn	-= Eve;

	return (fn);
}

double fn_Evib::df(double T)
{
	double dfn = 0.0;

	for (int s = 0; s <= NS-1; s++)
	{
		int ss = species.whereis[s];
		dfn += rhos[s]* (*species.data_entire)[ss].Cv_VEE(T);
	}

	return (dfn);
}



fn_Evib2::fn_Evib2()
{
	NS	= 0;
	Eve	= 0;
}

fn_Evib2::~fn_Evib2()
{

}

void fn_Evib2::create(SPECIES_DATA_BASIC &species_data, vector<double> &rho, int ns, double E)
{
	species	= species_data;
	rhos	= rho;
	NS		= ns;
	Eve		= E;
}

double fn_Evib2::f(double T)
{
	double fn	= 0.0;

	for (int s = 0; s <= NS-1; s++)
	{
		int ss = species.whereis[s];
		fn += rhos[s]* (*species.data_entire)[ss].e_VE(T);
	}

	fn	-= Eve;

	return (fn);
}

double fn_Evib2::df(double T)
{
	double dfn = 0.0;

	for (int s = 0; s <= NS-1; s++)
	{
		int ss = species.whereis[s];
		dfn += rhos[s]* (*species.data_entire)[ss].Cv_VE(T);
	}

	return (dfn);
}
