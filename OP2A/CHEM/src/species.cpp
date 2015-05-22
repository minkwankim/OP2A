/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: May 27, 2014
 *      			Author: Minkwan Kim
 *
 * species.cpp
 * 			-  
 *  
 */

#include "../include/species.hpp"
#include "../include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"


/*
 * Class NUM 4
 * 	- Species data
 */
Species::Species()
{

}

Species::~Species()
{

}

void Species::asg_data(string name_sp)
{
	basic.name	= name_sp;
	basic.read();

	thermo.assign_data(basic);
	transport.assign_data(basic);
}





double Species::Cv_vib(double T)
{
	int n;
	double Cv;
	double aux2, aux3;

	Cv	= 0.0;
	if (basic.type == MOLECULE)
	{
		for (n = 0;  n <= thermo.n_vib_lvl-1; n++)
		{
			aux2	= thermo.theta_vib[n]/T;
			aux3	= exp(aux2);
			Cv		+= thermo.R * aux2*aux2*aux3 / pow((aux3-1.0), 2.0);
		}
	}

	return (Cv);
}

double Species::Cv_ele(double T)
{
	int n;
	double Cv;
	double aux2, aux3, aux4;
	double component1, component2, component3, component4;

	component1	= 0.0;
	component2	= 0.0;
	component3	= 0.0;
	component4	= 0.0;

	Cv	= 0.0;
	if (basic.type != ELECTRON)
	{
		aux2		= thermo.theta_el[0]/ T;
		aux3		= exp(-aux2);

		component1	+= thermo.g[0] * aux3;
		component2 	+= thermo.g[0] * aux2/T * aux3;

		for (n = 1;  n <= thermo.n_elec_lvl; n++)
		{
			aux2		= thermo.theta_el[n]/ T;
			aux3		= exp(-aux2);

			component1	+= thermo.g[n] * aux3;
			component2 	+= thermo.g[n] * aux2/T * aux3;
			component3	+= thermo.g[n] * aux2*aux2 * aux3;
			component4  += thermo.g[n] * thermo.theta_el[n] * aux3;
		}

		Cv	= thermo.R * ((component3/component1) - (component4*component2)/(component1*component1));
	}

	return (Cv);
}

double Species::Cv_VE(double T, int NONEQ_VIB, int NONEQ_E)
{
	double Cv;

	switch(NONEQ_VIB)
	{
	case 0:
		Cv	= Cv_vib(T)	+ Cv_ele(T);
		break;
	case 1:
		Cv	= Cv_vib(T);
		break;
	case 2:
		Cv	= Cv_vib(T);
		break;
	case 3:
		if (basic.type	== ELECTRON)	Cv	= thermo.Cv[0];
		else							Cv	= Cv_vib(T)	+ Cv_ele(T);
		break;
	case 4:
		if (basic.type	== ELECTRON)	Cv	= thermo.Cv[0];
		else							Cv	= Cv_vib(T)	+ Cv_ele(T);
		break;
	}

	if (NONEQ_E != 0 && basic.type	== ELECTRON)	Cv = 0.0;

	return (Cv);
}


double Species::e_vib(double T)
{
	int n;
	double e;

	e	= 0.0;
	if (basic.type == MOLECULE)
	{
		for (n = 0;  n <= thermo.n_vib_lvl-1; n++)
		{
			e	+= thermo.R * thermo.theta_vib[n] / (exp(thermo.theta_vib[n]/T) - 1.0);
		}
	}

	return(e);
}


double Species::e_ele(double T)
{
	int n;
	double e_el;
	double aux2, aux3, aux4;

	e_el	= 0.0;
	if (basic.type != ELECTRON)
	{
		aux2	= 0.0;
		aux3	= 0.0;

		if (thermo.n_elec_lvl > 0)
		{
			for (n = 0;  n <= thermo.n_elec_lvl-1; n++)
			{
				aux4	= thermo.g[n]	* exp(-thermo.theta_el[n]/ T);
				aux2	+= aux4 * thermo.theta_el[n];
				aux3	+= aux4;
			}

			e_el	= thermo.R * aux2/aux3;
		}
	}

	return (e_el);
}

double Species::e_VE(double T, int NONEQ_VIB, int NONEQ_E)
{
	double e;

	switch(NONEQ_VIB)
	{
	case 0:
		e	= e_vib(T)	+ e_ele(T);
		break;
	case 1:
		e	= e_vib(T);
		break;
	case 2:
		e	= e_vib(T);
		break;
	case 3:
		if (basic.type	== ELECTRON)	e	= thermo.Cv[0] * T;
		else							e	= e_vib(T)	+ e_ele(T);
		break;
	case 4:
		if (basic.type	== ELECTRON)	e	= thermo.Cv[0] * T;
		else							e	= e_vib(T)	+ e_ele(T);
		break;
	}

	if (NONEQ_E	!= 0)	e = 0.0;

	return (e);
}

double Species::e_internal(double T, double Tr, double Tv, double Te, int NONEQ_ROT, int NONEQ_VIB, int NONEQ_E)
{
	double es;


	if (basic.type != ELECTRON)
	{
		if(NONEQ_VIB !=	0)
		{
			es	= thermo.Cv[0]*T	+ thermo.Cv[1]*Tr	+ e_VE(Tv, NONEQ_VIB, NONEQ_E)	+ basic.h0;
		}
		else
		{
			es	= thermo.Cv[0]*T	+ thermo.Cv[1]*Tr	+ basic.h0;
		}
	}
	else
	{
		es	= thermo.Cv[0] * Te;
	}

	return (es);
}

double Species::h(double T, double Tr, double Tv, double Te, int NONEQ_ROT, int NONEQ_VIB, int NONEQ_E)
{
	double hs;
	double es;

	es	= e_internal(T, Tr, Tv, Te, NONEQ_ROT, NONEQ_VIB, NONEQ_E);

	if (basic.type != ELECTRON)
	{
		hs	= thermo.R*T	+ es;
	}
	else
	{
		hs	= thermo.R*Te	+ es;
	}

	return (hs);
}


