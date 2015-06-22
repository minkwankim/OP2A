/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 19, 2015
 *      			Author: Minkwan Kim
 *
 * CalulateIC.cpp
 * 			-  
 *  
 */

#include "Common/include/CodeLocation.hpp"
#include "Common/include/Exception_NPExceed.hpp"
#include "Common/include/StringOps.hpp"
#include "CFD/include/VariableSamplesWIG.hpp"

#include "CFD/include/VariableConstants.hpp"
#include "CFD/include/VariableChange.hpp"
#include "CFD/include/VariableSet.hpp"
#include "Math/include/MathMisc.hpp"

#include "../include/OP2A_Application.hpp"



void ApplicationOP2A::CalculateIC()
{
	// 1. Allocate data of ICs
	IC_Q.resize(problem_setup.IC.NIC);
	IC_V.resize(problem_setup.IC.NIC);

	for (int n = 0; n <= problem_setup.IC.NIC-1; n++)
	{
		IC_Q(n)	= CFD::CFD_VariableSet::Q(species_set, grid.ND, problem_setup.NER, problem_setup.NEV, problem_setup.NEE, problem_setup.is_viscous, problem_setup.is_axisymmetric);
		IC_V(n)	= CFD::CFD_VariableSet::V(species_set, grid.ND, problem_setup.NER, problem_setup.NEV, problem_setup.NEE, problem_setup.is_viscous, problem_setup.is_axisymmetric);
	}


	// 2. Allocate V data
	int index;
	string var_name_temp;

#pragma omp parallel for num_threads(NT) private(var_name_temp)
	for (int n = 0; n <= problem_setup.IC.NIC-1; n++)
	{
		for (int s = 0; s <= species_set.NS-1; s++)
		{
			var_name_temp	= CFD::CFD_VariableSet::var_stringV(species_set.species[s].name, CFD::FluxCategory::Mass);
			IC_V.data[n](var_name_temp) =  problem_setup.IC.rho_s[n][s];
		}


		for (int k = 0; k <= grid.ND-1; k++)
		{
			var_name_temp	= CFD::CFD_VariableSet::var_stringV(k, CFD::FluxCategory::Momentum);
			IC_V.data[n](var_name_temp) =  problem_setup.IC.v[n][k];
		}


		var_name_temp	=  CFD::CFD_VariableSet::var_stringV(0, CFD::FluxCategory::Energy);
		IC_V.data[n](var_name_temp) =  problem_setup.IC.T[n];

		if (problem_setup.NER != 0)
		{
			var_name_temp	=  CFD::CFD_VariableSet::var_stringV(CHEM::EnergyMode::E_ROT, CFD::FluxCategory::Energy);
			IC_V.data[n](var_name_temp) =  problem_setup.IC.Tr[n];
		}


		if (problem_setup.NEV != 0)
		{
			var_name_temp	=  CFD::CFD_VariableSet::var_stringV(CHEM::EnergyMode::E_VE, CFD::FluxCategory::Energy);
			IC_V.data[n](var_name_temp) =  problem_setup.IC.Tv[n];
		}


		if (problem_setup.NEE != 0)
		{
			var_name_temp	=  CFD::CFD_VariableSet::var_stringV(CHEM::EnergyMode::E_ELE, CFD::FluxCategory::Energy);
			IC_V.data[n](var_name_temp) =  problem_setup.IC.Te[n];
		}
	}

	int nt = Math::fmin<int>(problem_setup.IC.NIC, NT);

#pragma omp parallel for num_threads(nt)
	for (int n = 0; n <= problem_setup.IC.NIC-1; n++)
	{
		CFD::VariableChange::V_to_Q(CFD_variabletype, CFD_NT, IC_V.data[n], species_set, grid.ND, IC_Q.data[n]);
	}
}
