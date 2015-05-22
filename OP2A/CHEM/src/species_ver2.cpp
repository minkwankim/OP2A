/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 15, 2015
 *      			Author: Minkwan Kim
 *
 * species_class.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_error_message.hpp"



// Constructor and Destructor
SPECIES_ver2::SPECIES_ver2()
{

}

SPECIES_ver2::~SPECIES_ver2()
{

}


/*
 * Internal Functions
 */

// F01: Assign Data
void SPECIES_ver2::asg_data(string name_sp)
{
	basic_data.name	= name_sp;
	basic_data.read();

	thermodynamic_properties.assign_data(basic_data);
	transport_properties.assign_data(basic_data);
}


// F02: Vibrational specific heat calculation
double SPECIES_ver2::Cv_vib(double T)
{
	int n;
	double Cv;
	double aux2, aux3;

	Cv	= 0.0;
	if (basic_data.type == MOLECULE)
	{
		for (n = 0;  n <= thermodynamic_properties.n_vib_lvl-1; n++)
		{
			aux2	= thermodynamic_properties.theta_vib[n]/T;
			aux3	= exp(aux2);
			Cv		+= thermodynamic_properties.R * aux2*aux2*aux3 / pow((aux3-1.0), 2.0);
		}
	}

	return (Cv);
}


// F03:: Electronic Specific heat calculation
double SPECIES_ver2::Cv_ele(double T)
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
	if (basic_data.type != ELECTRON)
	{
		aux2		= thermodynamic_properties.theta_el[0]/ T;
		aux3		= exp(-aux2);

		component1	+= thermodynamic_properties.g[0] * aux3;
		component2 	+= thermodynamic_properties.g[0] * aux2/T * aux3;

		for (n = 1;  n <= thermodynamic_properties.n_elec_lvl; n++)
		{
			aux2		= thermodynamic_properties.theta_el[n]/ T;
			aux3		= exp(-aux2);

			component1	+= thermodynamic_properties.g[n] * aux3;
			component2 	+= thermodynamic_properties.g[n] * aux2/T * aux3;
			component3	+= thermodynamic_properties.g[n] * aux2*aux2 * aux3;
			component4  += thermodynamic_properties.g[n] * thermodynamic_properties.theta_el[n] * aux3;
		}

		Cv	= thermodynamic_properties.R * ((component3/component1) - (component4*component2)/(component1*component1));
	}

	return (Cv);
}


// F03:: Vibrational-Electronic Specific heat calculation
double SPECIES_ver2::Cv_VE(double T)
{
	double Cv;
	Cv	= Cv_vib(T);//	+ Cv_ele(T);
	return (Cv);
}


// F04:: Vibrational-Electronic-Electron Specific heat calculation
double SPECIES_ver2::Cv_VEE(double T)
{
	double Cv;

	if (basic_data.type != ELECTRON)	Cv	= Cv_vib(T)	+ Cv_ele(T);
	else								Cv	= thermodynamic_properties.Cv[TRA];
	return (Cv);
}



// F05:: Vibrational specific Energy
double SPECIES_ver2::e_vib(double T)
{
	int n;
	double e;

	e	= 0.0;
	if (basic_data.type == MOLECULE)
	{
		for (n = 0;  n <= thermodynamic_properties.n_vib_lvl-1; n++)
		{
			e	+= thermodynamic_properties.R * thermodynamic_properties.theta_vib[n] / (exp(thermodynamic_properties.theta_vib[n]/T) - 1.0);
		}
	}

	return(e);
}


// F06:: Electronic specific Energy
double SPECIES_ver2::e_ele(double T)
{
	int n;
	double e_el;
	double aux2, aux3, aux4;

	e_el	= 0.0;
	if (basic_data.type != ELECTRON)
	{
		aux2	= 0.0;
		aux3	= 0.0;

		if (thermodynamic_properties.n_elec_lvl > 0)
		{
			for (n = 0;  n <= thermodynamic_properties.n_elec_lvl-1; n++)
			{
				aux4	= thermodynamic_properties.g[n]	* exp(-thermodynamic_properties.theta_el[n]/ T);
				aux2	+= aux4 * thermodynamic_properties.theta_el[n];
				aux3	+= aux4;
			}

			e_el	= thermodynamic_properties.R * aux2/aux3;
		}
	}

	return (e_el);
}


// F07:: Vibrational-electronic Specific energy
double SPECIES_ver2::e_VE(double T)
{
	double e_ve;

	e_ve	= e_vib(T);//	+ e_ele(T);
	return (e_ve);
}


// F08:: Vibrational-electronic-electron Specific energy
double SPECIES_ver2::e_VEE(double T)
{
	double e_vee;

	if (basic_data.type != ELECTRON)	e_vee	= e_VE(T);
	else								e_vee	= 0.0; //Cv_VEE(T) * T;
	return (e_vee);
}



// F09: Internal energy per unit mass
double SPECIES_ver2::e_internal(double T, double Tr, double Tv, double Tele)
{

	double es;

	es	= thermodynamic_properties.Cv[0]*T	+ thermodynamic_properties.Cv[1]*Tr	+ basic_data.h0;
	es	+= e_vib(Tv) + e_ele(Tele);

	return (es);
}



// F10: species enthalpy
double SPECIES_ver2::h(double T, double Tr, double Tv, double Tele)
{
	double hs;
	double es;

	hs	= e_internal(T, Tr, Tv, Tele)	+ thermodynamic_properties.R*T;
	return (hs);
}



double SPECIES_ver2::tau_e_vib(double Te, double ne)
{
	int i;
	double v0, v, a, b, c, x;
	double tau_ev, k_eV;
	double temp1, temp2, dv;


	temp1 = 0.5 * ne * pow((1.0 - exp(-thermodynamic_properties.theta_vib[0]/Te)), 2.0);
	temp2 = 0.0;

	v0 = 0.0;
	for (i = 0; i <= thermodynamic_properties.n_vib_ele-1; i++)
	{
		v = thermodynamic_properties.kev[i][0];
		a = thermodynamic_properties.kev[i][1];
		b = thermodynamic_properties.kev[i][2];
		c = thermodynamic_properties.kev[i][3];

		x = Te * K_TO_eV;
		k_eV = 1.0e-15 * a*pow(x, -1.5) * exp(b/x + c);

		dv = fabs(v - v0);
		v0 = v;

		temp2 += k_eV*pow(v,2.0)*dv;
	}

	tau_ev = temp1*temp2;
	if (error_check_double_pos(tau_ev) == true)
	{
		Error_message_type	error_message;
		error_message.module_name	="CHEMISTRY MODULE";
		error_message.location_primary_name	= "N/A";
		error_message.location_secondary_name	= "N/A";
		error_message.message	= "Problem in cal_tau_e_vib.c: It gives negative relaxation time";
		error_message.print_message();
	}

	tau_ev = 1.0/tau_ev;
	return(tau_ev);
}
