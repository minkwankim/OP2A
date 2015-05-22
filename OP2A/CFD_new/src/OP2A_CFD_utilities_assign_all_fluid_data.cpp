/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 17, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_utilities_assign_all_fluid_data.cpp
 * 			-  
 *  
 */




#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include "../../CHEM/include/OP2A_chemistry.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
using namespace std;


void OP2A_CFD_assign_rho_and_T_for_all_fluid(int NCM, vector<SOL_CFD> &Solutions, vector<SPECIES_DATA_BASIC>	&species,
											vector<vector <double *> > &rho_s_ALL, vector<vector<vector <double *> > > &Ts_ALL, int num_fluid)
{
	for (int f = 0; f <= num_fluid-1; f++)
	{
		int ne	= Solutions[f].setup.NS + Solutions[f].setup.ND;

		for (int s = 0; s <= Solutions[f].setup.NS-1; s++)
		{
			int s_ptr	= species[f].whereis[s];
			for (int c = 0; c <= NCM-1; c++)	rho_s_ALL[c][s_ptr]	= &Solutions[f].Vc[c][s];


			switch ((*species[f].data_entire)[s_ptr].basic_data.type)
			{
			case MOLECULE:
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][TRA]	= &Solutions[f].Vc[c][ne];
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][ROT]	= &Solutions[f].Vc[c][ne+Solutions[f].setup.ID_T[ROT]];
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][VIB]	= &Solutions[f].Vc[c][ne+Solutions[f].setup.ID_T[VIB]];
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][ELE]	= new double[1];

				break;

			case ATOM:
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][TRA]	= &Solutions[f].Vc[c][ne];
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][ROT]	= new double[1];
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][VIB]	= new double[1];
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][ELE]	= new double[1];
				break;

			case ELECTRON:
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][TRA]	= &Solutions[f].Vc[c][ne];
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][ROT]	= new double[1];
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][VIB]	= new double[1];
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][ELE]	= &Solutions[f].Vc[c][ne+Solutions[f].setup.ID_T[ELE]];
				break;
			}
		}
	}
}


void OP2A_CFD_assign_rho_and_T_for_all_fluid_ver2(int NCM, vector<SOL_CFD> &Solutions, vector<SPECIES_DATA_BASIC>	&species,
												vector<vector <double> > &rho_s_ALL, vector<vector<vector <double> > > &Ts_ALL, int num_fluid)
{
	for (int f = 0; f <= num_fluid-1; f++)
	{
		int ne	= Solutions[f].setup.NS + Solutions[f].setup.ND;

		for (int s = 0; s <= Solutions[f].setup.NS-1; s++)
		{
			int s_ptr	= species[f].whereis[s];
			for (int c = 0; c <= NCM-1; c++)	rho_s_ALL[c][s_ptr]	= Solutions[f].Vc[c][s];


			switch ((*species[f].data_entire)[s_ptr].basic_data.type)
			{
			case MOLECULE:
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][TRA]	= Solutions[f].Vc[c][ne];
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][ROT]	= Solutions[f].Vc[c][ne+Solutions[f].setup.ID_T[ROT]];
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][VIB]	= Solutions[f].Vc[c][ne+Solutions[f].setup.ID_T[VIB]];
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][ELE]	= 0.0;

				break;

			case ATOM:
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][TRA]	= Solutions[f].Vc[c][ne];
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][ROT]	= 0.0;
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][VIB]	= 0.0;
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][ELE]	= 0.0;
				break;

			case ELECTRON:
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][TRA]	= Solutions[f].Vc[c][ne];
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][ROT]	= 0.0;
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][VIB]	= 0.0;
				for (int c = 0; c <= NCM-1; c++)	Ts_ALL[c][s_ptr][ELE]	= Solutions[f].Vc[c][ne+Solutions[f].setup.ID_T[ELE]];
				break;
			}
		}
	}
}

