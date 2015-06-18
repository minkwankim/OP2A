/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 15, 2015
 *      			Author: Minkwan Kim
 *
 * ApplicationInitializeData.cpp
 * 			-  
 *  
 */


#include "Common/include/CodeLocation.hpp"
#include "Common/include/Exception_NPExceed.hpp"
#include "Common/include/StringOps.hpp"
#include "CFD/include/VariableSamplesWIG.hpp"
#include "CFD/include/VariableConstants.hpp"



//#include "../include/ApplicationConstants.hpp"
#include "../include/OP2A_Application.hpp"






void ApplicationOP2A::InitializeData(unsigned int num_ic)
{
	int index;
	string var_name_temp;

	for (int s = 0; s <= species_set.NS-1; s++)
	{
		var_name_temp	= CFD::CFD_VariableSet::var_stringV(species_set.species[s].name, CFD::FluxCategory::Mass);

#		pragma omp parallel for num_threads(NT)
		for (int c = 0; c <= grid.NCM; c++)
		{
			grid.cells[c].data1D(NAME_V, var_name_temp)	= problem_setup.IC.rho_s[num_ic][s];
		}
	}


	for (int k = 0; k <= grid.ND-1; k++)
	{
		var_name_temp	= CFD::CFD_VariableSet::var_stringV(k, CFD::FluxCategory::Momentum);

#		pragma omp parallel for num_threads(NT)
		for (int c = 0; c <= grid.NCM; c++)
		{
			grid.cells[c].data1D(NAME_V, var_name_temp)	= problem_setup.IC.v[num_ic][k];
		}
	}


	var_name_temp	=  CFD::CFD_VariableSet::var_stringV(0, CFD::FluxCategory::Energy);
#	pragma omp parallel for num_threads(NT)
	for (int c = 0; c <= grid.NCM; c++)
	{
		grid.cells[c].data1D(NAME_V, var_name_temp)	= problem_setup.IC.T[num_ic];
	}


	if (problem_setup.NER != 0)
	{
		var_name_temp	=  CFD::CFD_VariableSet::var_stringV(CHEM::EnergyMode::E_ROT, CFD::FluxCategory::Energy);

#		pragma omp parallel for num_threads(NT)
		for (int c = 0; c <= grid.NCM; c++)
		{
			grid.cells[c].data1D(NAME_V, var_name_temp)	= problem_setup.IC.Tr[num_ic];
		}
	}


	if (problem_setup.NEV != 0)
	{
		var_name_temp	=  CFD::CFD_VariableSet::var_stringV(CHEM::EnergyMode::E_VIB, CFD::FluxCategory::Energy);

#		pragma omp parallel for num_threads(NT)
		for (int c = 0; c <= grid.NCM; c++)
		{
			grid.cells[c].data1D(NAME_V, var_name_temp)	= problem_setup.IC.Tv[num_ic];
		}
	}


	if (problem_setup.NEE != 0)
	{
		var_name_temp	=  CFD::CFD_VariableSet::var_stringV(CHEM::EnergyMode::E_ELE, CFD::FluxCategory::Energy);

#		pragma omp parallel for num_threads(NT)
		for (int c = 0; c <= grid.NCM; c++)
		{
			grid.cells[c].data1D(NAME_V, var_name_temp)	= problem_setup.IC.Te[num_ic];
		}
	}
}

