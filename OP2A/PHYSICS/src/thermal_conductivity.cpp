/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * thermal_conductivity.cpp
 * 			-  
 *  
 */


#include "../include/constant_text.hpp"
#include "../include/transport_properties.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"


double thermal_conductivity_Eucken(double mu_s, double T, Species* species_single, int MODE)
{
	double kappa;
	double Cv;
	switch (MODE)
	{
	case TRA:
		Cv		= species_single->thermo.Cv[0];
		kappa	= 2.5*mu_s*Cv;
		break;

	case ROT:
		Cv		= species_single->thermo.Cv[1];
		kappa	= mu_s*Cv;
		break;

	case VIB:
		Cv		= species_single->Cv_vib(T);
		kappa	= mu_s*Cv;
		break;

	case ELE:
		Cv		= species_single->thermo.Cv[0];
		kappa	= mu_s*Cv;
		break;

	case VE:
		Cv		= species_single->Cv_vib(T) + species_single->Cv_ele(T);
		kappa	= mu_s*Cv;
		break;

	case TR:
		kappa	= mu_s * (2.5*species_single->thermo.Cv[0] + species_single->thermo.Cv[1]);
		break;
	}

	return (kappa);
}



double calculate_thermal_conductivity_Eucken(double mu_s, double T, SPECIES species, int MODE)
{
	double kappa;
	double Cv;

	switch (MODE)
	{
	case 0:	// Translation
		Cv	= species.thermodynamic_properties.Cv[0];
		Cv	= 2.5 *Cv;
		break;
	case 1:	// Rotation;
		Cv	= species.thermodynamic_properties.Cv[0];
		break;
	case 2:	// Vibration
		Cv	= species.Cv_VE(T);
	}

	kappa	= Cv * mu_s;

	return (kappa);
}


void calculate_thermal_conductivity(double mu_s, double Tv, SPECIES species, int method, vector<double> &kappa)
{
	double kappa_tra	= 0.0;
	double kappa_rot	= 0.0;
	double kappa_vib	= 0.0;

	switch (method)
	{
	case 0:
		kappa_tra	= calculate_thermal_conductivity_Eucken(mu_s, Tv, species, 0);

		if (species.basic_data.type == MOLECULE)
		{
			kappa_rot	= calculate_thermal_conductivity_Eucken(mu_s, Tv, species, 1);
			kappa_vib	= calculate_thermal_conductivity_Eucken(mu_s, Tv, species, 2);
		}

		break;
	}


	kappa[TRA]	= kappa_tra;
	kappa[ROT]	= kappa_rot;
	kappa[VIB]	= kappa_vib;
}
