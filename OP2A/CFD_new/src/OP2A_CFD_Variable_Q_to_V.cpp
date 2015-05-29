/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Variable_Q_to_V.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_CFD_Variable_change.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../../MATRIX/include/Newton_solver.hpp"


void CFD_Q_to_V(vector<double> Q, vector<double> &V, CFD_variable_setup_ver2 setup, CFD_mixture_data variable_data, SPECIES_DATA_BASIC &species)
{
	double epsilon	= 1.0e-3;


	if (variable_data.rho != 0.0)
	{
#pragma omp parallel for
		for (int s = 0; s <= setup.NS-1; s++)	V[s]			= Q[s];

#pragma omp parallel for
		for (int k = 0; k <= setup.ND-1; k++)	V[setup.NS+k]	= Q[setup.NS+k] / variable_data.rho;

		int ne		= setup.NS+setup.ND;
		double E 		= Q[ne];
		double rho_Cv	= variable_data.rho * variable_data.Cv_bar;

		for (int s = 0; s <= setup.NS-1; s++)	E -= Q[s]*species.h0[s];
		for (int k = 0; k <= setup.ND-1; k++)	E -= 0.5*variable_data.rho*V[setup.NS+k]*V[setup.NS+k];
		for (int k = 1; k <= setup.NE-1; k++)	E -= Q[ne+k];

		if (E <= 0)
		{
			Error_message_type error;
			error.module_name	= "CFD module: Q_to_V:";
			error.location_primary_name	= "N/A";
			error.location_secondary_name = "N/A";
			error.message = "Negative temperature.";
			error.print_message();
		}
		else
		{
			V[ne]	= E / rho_Cv;
		}


		fn_Evib	fn;
		vector<double>	rhos (setup.NS, 0.0);
		Newton_solver<fn_Evib>	solve_Tve;

		if (setup.energy_flag == 1)
		{
			V[ne+1]	= Q[ne+1] / (species.Cv_t[setup.NS-1]*Q[setup.NS-1]);
		}
		else if (setup.energy_flag == 2)
		{
			vector<double>	rhos (setup.NS, 0.0);
			for (int s = 0; s <= setup.NS-1; s++)	rhos[s]	=  Q[s];
			fn.create(species, rhos, setup.NS, Q[ne+1]);

			Newton_solver<fn_Evib>	solve_Tve;
			solve_Tve.Iter_MAX					= 10000;
			solve_Tve.epsilon					= epsilon;
			solve_Tve.X_min						= 10.0;
			solve_Tve.X_max						= 200000.0;
			V[ne+1]	=	solve_Tve.solve(fn, 1);
		}
		else if (setup.energy_flag == 3)
		{
			for (int s = 0; s <= setup.NS-1; s++)	rhos[s]	=  Q[s];
			fn.create(species, rhos, setup.NS, Q[ne+1]);

			solve_Tve.Iter_MAX					= 10000;
			solve_Tve.epsilon					= epsilon;
			solve_Tve.X_min						= 10.0;
			solve_Tve.X_max						= 200000.0;

			V[ne+1]	= solve_Tve.solve(fn, 1);
			V[ne+2]	= Q[ne+2] / (species.Cv_t[setup.NS-1]*Q[setup.NS-1]);
		}
		else if (setup.energy_flag == 4)
		{
			V[ne+1]	= Q[ne+1] / (variable_data.rho*variable_data.Cv_bar_rot);
		}
		else if (setup.energy_flag == 5)
		{
			V[ne+1]	= Q[ne+1] / (variable_data.rho*variable_data.Cv_bar_rot);
			V[ne+2]	= Q[ne+2] / (species.Cv_t[setup.NS-1]*Q[setup.NS-1]);
		}
		else if (setup.energy_flag == 6)
		{
			V[ne+1]	= Q[ne+1] / (variable_data.rho*variable_data.Cv_bar_rot);

			for (int s = 0; s <= setup.NS-1; s++)	rhos[s]	=  Q[s];
			fn.create(species, rhos, setup.NS, Q[ne+2]);

			solve_Tve.Iter_MAX					= 10000;
			solve_Tve.epsilon					= epsilon;
			solve_Tve.X_min						= 10.0;
			solve_Tve.X_max						= 200000.0;

			V[ne+2]	= solve_Tve.solve(fn, 1);
		}
		else if (setup.energy_flag == 7)
		{
			V[ne+1]	= Q[ne+1] / (variable_data.rho*variable_data.Cv_bar_rot);

			for (int s = 0; s <= setup.NS-1; s++)	rhos[s]	=  Q[s];
			fn.create(species, rhos, setup.NS, Q[ne+2]);

			solve_Tve.Iter_MAX					= 10000;
			solve_Tve.epsilon					= epsilon;
			solve_Tve.X_min						= 10.0;
			solve_Tve.X_max						= 200000.0;

			V[ne+2]	= solve_Tve.solve(fn, 1);
			V[ne+3]	= Q[ne+3] / (species.Cv_t[setup.NS-1]*Q[setup.NS-1]);
		}
	}
	else
	{
#pragma omp parallel for
		for (int i = 0; i <= setup.VAR-1; i++)	V[i]	= 0.0;
	}
}
