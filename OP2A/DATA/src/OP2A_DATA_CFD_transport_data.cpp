/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_DATA_CFD_transport_data.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_DATA_CFD_transport_data.hpp"
#include "../../PHYSICS/include/transport_properties.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"


CFD_transport_properties_data::CFD_transport_properties_data()
{
	viscosity_coefficient_mixture	= 0.0;
}



CFD_transport_properties_data::~CFD_transport_properties_data()
{

}

void CFD_transport_properties_data::allocate_size(int ns, int ne)
{
	viscosity_coefficient_mixture	= 0.0;

	diffusion_coefficient.resize(ns);
	viscosity_coefficient.resize(ns);

	thermal_conductivity_mixture.resize(ne);
	thermal_conductivity = vector_2D<double>(ns, ne, 0.0);
}



void CFD_transport_properties_data::mean_value(CFD_transport_properties_data &data1, CFD_transport_properties_data &data2, int NS, int NE)
{
	for (int s = 0; s <= NS-1; s++)	diffusion_coefficient[s]	= 0.5 * (data1.diffusion_coefficient[s] + data2.diffusion_coefficient[s]);
	for (int s = 0; s <= NS-1; s++)	viscosity_coefficient[s]	= 0.5 * (data1.viscosity_coefficient[s] + data2.viscosity_coefficient[s]);

	viscosity_coefficient_mixture	= 0.5 * (data1.viscosity_coefficient_mixture + data2.viscosity_coefficient_mixture);


#pragma omp parallel
	for (int s = 0; s <= NS-1; s++)
	{
		for (int k = 0; k <= NE-1; k++)
		{
			thermal_conductivity[s][k]	= 0.5 * (data1.thermal_conductivity[s][k] + data2.thermal_conductivity[s][k]);
		}
	}

#pragma omp parallel for
	for (int k = 0; k <= NE-1; k++)
	{
		thermal_conductivity_mixture[k]	= 0.5 * (data1.thermal_conductivity_mixture[k] + data2.thermal_conductivity_mixture[k]);
	}
}


void CFD_transport_properties_data:: calculate_transport_properties(CFD_variable_setup_ver2 &setup, vector<double> &V, vector<double> &Xs, SPECIES_DATA_BASIC &species, CFD_transport_properties_setup_ver2 viscosity_setup)
{
	double T;
	double Tv;

	T 	= V[setup.NS+setup.NE];
	Tv	= V[setup.NS+setup.NE+setup.ID_T[VIB]];


	// 1. Calculate viscosity
#pragma omp parallel for
	for (int s = 0; s <= setup.NS-1; s++)
	{
		double mu;
		int s_ptr	= species.whereis[s];

		mu	= viscosity_mu(T, viscosity_setup.viscosity_model, (*species.data_entire)[s_ptr]);
		viscosity_coefficient[s]	= mu;
	}


		// 2. Conductivity
	vector<double>	kappa(3, 0.0);

#pragma omp parallel for num_threads(species.NS)
	for (int s = 0; s <= species.NS-1; s++)
	{
		int s_ptr	= species.whereis[s];
		calculate_thermal_conductivity(viscosity_coefficient[s], Tv, (*species.data_entire)[s_ptr], viscosity_setup.conductivity_model, kappa);

		thermal_conductivity[s][0] = kappa[0];


		switch(setup.energy_flag)
		{
		case 0:
			thermal_conductivity[s][0] += kappa[ROT];
			thermal_conductivity[s][0] += kappa[VIB];
			break;
		case 1:
			thermal_conductivity[s][0] += kappa[ROT];
			thermal_conductivity[s][0] += kappa[VIB];
			break;
		case 2:
			thermal_conductivity[s][0] += kappa[ROT];
			thermal_conductivity[s][1] = kappa[VIB];
			break;
		case 3:
			thermal_conductivity[s][0] += kappa[ROT];
			thermal_conductivity[s][1] = kappa[VIB];
			break;
		case 4:
			thermal_conductivity[s][1] = kappa[ROT];
			thermal_conductivity[s][0] += kappa[VIB];
			break;
		case 5:
			thermal_conductivity[s][1] = kappa[ROT];
			thermal_conductivity[s][0] += kappa[VIB];
			break;
		case 6:
			thermal_conductivity[s][1] = kappa[ROT];
			thermal_conductivity[s][1] = kappa[VIB];
			break;
		case 7:
			thermal_conductivity[s][1] = kappa[ROT];
			thermal_conductivity[s][1] = kappa[VIB];
			break;
		}
	}



	// 3. Mixture properties
	switch (viscosity_setup.mixing_rule)
	{
	case 0:
		vector<double> phi_s(species.NS, 0.0);
		for (int s = 0; s <= species.NS-1; s++)
		{
			phi_s[s]	= 0.0;
			for (int r = 0; r <= setup.NS-1; r++)
			{
				if ((*species.data_entire)[species.whereis[r]].basic_data.type != ELECTRON)
				{
					double temp1	= 1.0 + sqrt(viscosity_coefficient[s]/viscosity_coefficient[r]) * pow(species.M[r]/species.M[s], 0.25);
					double temp2	= 8.0 * ( 1.0 + species.M[s]/species.M[r]);
					phi_s[s]		+= Xs[r] * pow(temp1, 2.0) / sqrt(temp2);
				}
			}

			if(phi_s[s]	!= phi_s[s]	||  phi_s[s] < 0.0 || phi_s[s] == numeric_limits<double>::infinity())
			{
				Error_message_type	error_message;
				error_message.module_name	="Physics_properties: Wilke mixing rule";
				error_message.location_primary_name	= "N/A";
				error_message.location_secondary_name	= "N/A";
				error_message.message	= " Problem in the calculation of phi!";
				error_message.print_message();
			}
		}

		// 2. Mixture viscosity
		viscosity_coefficient_mixture  = 0.0;
		for (int s = 0; s <= setup.NS-1; s++)	viscosity_coefficient_mixture += Xs[s] * viscosity_coefficient[s] / phi_s[s];

		for (int k = 0; k <= setup.NE-1; k++)
		{
			thermal_conductivity_mixture[k]	= 0.0;
			for (int s = 0; s <= species.NS-1; s++)	thermal_conductivity_mixture[k] += Xs[s] * thermal_conductivity[s][k] / phi_s[s];
		}
		break;
	}
}



