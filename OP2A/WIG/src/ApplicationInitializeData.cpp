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



#include "../include/ApplicationConstants.hpp"
#include "../include/OP2A_Application.hpp"






void ApplicationOP2A::InitializeData(unsigned int num_ic)
{
	int index;
	string var_name_temp;

	for (int s = 0; s <= species_set.NS-1; s++)
	{
		var_name_temp	= VAR_RHO_PRE + species_set.species[s].name + VAR_RHO_POST;
		index			= s;
#pragma omp parallel for num_threads(NT)
		for (int c = 0; c <= grid.NCM; c++)
		{
			grid.cells[c].data1D(VAR_VECTOR_V, var_name_temp)	= problem_setup.IC.rho_s[num_ic][index];
		}
	}


	var_name_temp	= VAR_U;
	index			= 0;
#pragma omp parallel for num_threads(NT)
	for (int c = 0; c <= grid.NCM; c++)
	{
		grid.cells[c].data1D(VAR_VECTOR_V, var_name_temp)	= problem_setup.IC.v[num_ic][index];
	}


	var_name_temp	= VAR_V;
	index			= 1;
#pragma omp parallel for num_threads(NT)
	for (int c = 0; c <= grid.NCM; c++)
	{
		grid.cells[c].data1D(VAR_VECTOR_V, var_name_temp)	= problem_setup.IC.v[num_ic][index];
	}


	if (grid.ND == 3)
	{
		var_name_temp	= VAR_W;
		index			= 2;
#pragma omp parallel for num_threads(NT)
		for (int c = 0; c <= grid.NCM; c++)
		{
			grid.cells[c].data1D(VAR_VECTOR_V, var_name_temp)	= problem_setup.IC.v[num_ic][index];
		}
	}


	var_name_temp	= VAR_T;
#pragma omp parallel for num_threads(NT)
	for (int c = 0; c <= grid.NCM; c++)
	{
		grid.cells[c].data1D(VAR_VECTOR_V, var_name_temp)	= problem_setup.IC.T[num_ic];
	}


	if (problem_setup.NER != 0)
	{
		var_name_temp	= VAR_T_ROT;
#pragma omp parallel for num_threads(NT)
		for (int c = 0; c <= grid.NCM; c++)
		{
			grid.cells[c].data1D(VAR_VECTOR_V, var_name_temp)	= problem_setup.IC.Tr[num_ic];
		}
	}


	if (problem_setup.NEV != 0)
	{
		var_name_temp	= VAR_T_VIB;
#pragma omp parallel for num_threads(NT)
		for (int c = 0; c <= grid.NCM; c++)
		{
			grid.cells[c].data1D(VAR_VECTOR_V, var_name_temp)	= problem_setup.IC.Tv[num_ic];
		}
	}


	if (problem_setup.NEE != 0)
	{
		var_name_temp	= VAR_T_E;
#pragma omp parallel for num_threads(NT)
		for (int c = 0; c <= grid.NCM; c++)
		{
			grid.cells[c].data1D(VAR_VECTOR_V, var_name_temp)	= problem_setup.IC.Te[num_ic];
		}
	}





}
