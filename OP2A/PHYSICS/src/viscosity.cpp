/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * viscosity.cpp
 * 			-  
 *  
 */



#include "../include/constant_text.hpp"
#include "../include/transport_properties.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"



double viscosity(double T, int model, Species* species_single)
{
	double mu;
	double temp1, temp2, temp3;

	if (species_single->basic.type == ELECTRON)
	{
		mu	= 0.0;
	}
	else
	{
		if (model == VHS_MODEL)
		{
			temp1	= sqrt(PI*C_BOLTZMANN_SI*species_single->basic.m*species_single->thermo.T_ref);
			temp2	= 2.0*PI*pow(species_single->thermo.d_ref, 2.0) * (5.0-2.0*species_single->thermo.omega) * (7.0-2.0*species_single->thermo.omega);
			mu = 15.0 * temp1 / temp2 * pow(T/species_single->thermo.T_ref, species_single->thermo.omega);
		}
		else if (model	== BLOTTNER_MODEL)
		{
			double lnT = log(T);
			temp1	= species_single->transport.Blottner[0]*lnT + species_single->transport.Blottner[1];
			mu 		= 0.1*exp(temp1*lnT + species_single->transport.Blottner[2]);
		}
		else if (model	== 	SUTHERLAND_MODEL)
		{
			temp1 	= species_single->transport.Sutherland[1] + species_single->transport.Sutherland[2];
			temp2 	= T + species_single->transport.Sutherland[2];
			temp3 	= T / species_single->transport.Sutherland[1];

			mu		= species_single->transport.Sutherland[0] * (temp1/temp2) * pow(temp3, 1.5);
		}
		else if (model	== KINETIC_MODEL)
		{
			temp1	= 1.147 * pow(T/species_single->transport.T_eps, -0.145) + pow(T/species_single->transport.T_eps + 0.5, -2.0);
			mu		= (2.68E-6 * sqrt(species_single->basic.M*T)) / (pow(species_single->transport.sigma, 2.0)*temp1);
		}
	}

	return (mu);
}


double viscosity_collision(double T, double Te, double *gamms_s, Species *species_all, int NS)
{
	int		s, r;
	double	mu_mix	= 0.0;
	double	collision_term_sr;


	double temp1;
	temp1	= 0.0;


	for (s	= 0; s	<= NS-1; s++)
	{
		double aux1;

		if (species_all[s].basic.type != ELECTRON)
		{
			aux1	= 0.0;
			for (r = 0; r <= NS-1; r++)
			{
				if (species_all[r].basic.type != ELECTRON)
				{

					//collision_term_sr	= Delta_sr_2(T, species_all[s].M, species_all[s].M, double pi_omega_sr_22)
					//aux1	+= gamma_s[r] * Delta_sr_2
				}
			}
		}
	}


	return(mu_mix);
}







double viscosity_mu_ver1(double T, int model, SPECIES &species)
{
	double mu;
	double temp1, temp2, temp3;

	if (species.basic_data.type == ELECTRON)
	{
		mu	= 0.0;
	}
	else
	{
		if (model == VHS_MODEL)
		{
			temp1	= sqrt(PI*C_BOLTZMANN_SI*species.basic_data.m*species.thermodynamic_properties.T_ref);
			temp2	= 2.0*PI*pow(species.thermodynamic_properties.d_ref, 2.0) * (5.0-2.0*species.thermodynamic_properties.omega) * (7.0-2.0*species.thermodynamic_properties.omega);
			mu = 15.0 * temp1 / temp2 * pow(T/species.thermodynamic_properties.T_ref, species.thermodynamic_properties.omega);
		}
		else if (model	== BLOTTNER_MODEL)
		{
			double lnT = log(T);
			temp1	= species.transport_properties.Blottner[0]*lnT + species.transport_properties.Blottner[1];
			mu 		= 0.1*exp(temp1*lnT + species.transport_properties.Blottner[2]);
		}
		else if (model	== 	SUTHERLAND_MODEL)
		{
			temp1 	= species.transport_properties.Sutherland[1] + species.transport_properties.Sutherland[2];
			temp2 	= T + species.transport_properties.Sutherland[2];
			temp3 	= T / species.transport_properties.Sutherland[1];

			mu		= species.transport_properties.Sutherland[0] * (temp1/temp2) * pow(temp3, 1.5);
		}
		else if (model	== KINETIC_MODEL)
		{
			temp1	= 1.147 * pow(T/species.transport_properties.T_eps, -0.145) + pow(T/species.transport_properties.T_eps + 0.5, -2.0);
			mu		= (2.68E-6 * sqrt(species.basic_data.M*T)) / (pow(species.transport_properties.sigma, 2.0)*temp1);
		}
	}

	if (mu != mu	||  mu == numeric_limits<double>::infinity())
	{
		Error_message_type	error_message;
		error_message.module_name	="Physics_properties: Viscosity";
		error_message.location_primary_name	= "N/A";
		error_message.location_secondary_name	= "N/A";
		error_message.message	= " Problem in viscosity coefficient!";
		error_message.print_message();
	}

	return (mu);
}












