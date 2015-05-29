/*
 * plasma_physics_basic.cpp
 *
 *  Created on: Jan 16, 2014
 *      Author: minkwan
 */

#include "../include/plasma_physics_basic_parameter.hpp"


// F-01: Debye Length
//
double Debye_length_MKS(double ne, double Te)
{
	double debye_m;

	//debye_m = (EPS0_SI * C_BOLTZMANN_SI * Te) / (ne * ELECTRON_CHARGE_SI *ELECTRON_CHARGE_SI);
	debye_m = sqrt(debye_m);
	return (debye_m);
}

double Debye_length_Gaussian(double ne, double Te)
{
	double debye_cm;

	//debye_cm = (C_BOLTZMANN_CGS * Te) / (4*PI * ne * ELECTRON_CHARGE_CGS *ELECTRON_CHARGE_CGS);
	debye_cm = sqrt(debye_cm);
	return (debye_cm);
}

double Debye_length_eV(double ne_mks, double Te_eV)
{
	double debye_m;

	debye_m = (EPS0_SI * Te_eV) / (ne_mks * ELECTRON_CHARGE_SI);
	debye_m = sqrt(debye_m);
	return (debye_m);
}

// F-02: Thermal velocity
//
double V_thermal_MKS(double T, double m)
{
	double V_th;
	//V_th	= (2.0*C_BOLTZMANN_SI*T) / m;
	V_th	= sqrt(V_th);
	return (V_th);
}


// F-03: Plasma Frequency
double Omega_p_MKS(double ne)
{
	double omega_p;

	//omega_p	= (ne * ELECTRON_CHARGE_SI * ELECTRON_CHARGE_SI) / (Me * EPS0_SI);
	omega_p	= sqrt(omega_p);
	return (omega_p);
}

double Omega_p_CGS(double ne)
{
	double omega_p;

	//omega_p	= (4.0*PI*ne*ELECTRON_CHARGE_CGS*ELECTRON_CHARGE_CGS) / Me_cgs;
	omega_p	= sqrt(omega_p);
	return (omega_p);
}
