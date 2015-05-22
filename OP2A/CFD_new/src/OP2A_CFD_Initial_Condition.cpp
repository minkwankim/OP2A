/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Initial_Condition.cpp
 * 			-  
 *  
 */


#include <vector>
#include "../../GRID/include/OP2A_grid.hpp"
#include "../../DATA/include/OP2A_data.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../include/OP2A_CFD_Variable_change.hpp"


using namespace std;


void Calculate_IC_ver2(CFD_variable_setup_ver2	variable_setup,	vector<double> &rho_s, vector<double> &u, vector<double> &T,  vector<double> &Q, vector<double> &V, SPECIES_DATA_BASIC &species)
{
#pragma omp parallel for
	for (int s = 0; s <= variable_setup.NS-1; s++)	V[s]	= rho_s[s];

#pragma omp parallel for
	for (int k = 0; k <= variable_setup.ND-1; k++)	V[variable_setup.NS+k]	= u[k];

	int ne	= variable_setup.NS + variable_setup.ND;
#pragma omp parallel for
	for (int m = 0; m <= variable_setup.NE-1; m++)	V[ne+m]	= T[m];

	CFD_mixture_data variable_data;
	variable_data.calculate_data(variable_setup.NS, variable_setup.ID_T[ROT], V, species.Cv_t, species.Cv_r, species.R, species.M);
	CFD_V_to_Q(V, Q, variable_setup, variable_data, species);
}


void Calculate_IC(CFD_variable_setup_ver2	variable_setup,
				vector< vector < double > > &rho_s,
				vector< vector < double > > &u,
				vector <double> &T, vector <double> &Tr, vector <double> &Tv, vector <double> &Te,
				vector< vector<double> > &Q, vector< vector<double> > &V, SPECIES_DATA_BASIC &species)
{
	int N_IC	= T.size();
	Q	= vector_2D(N_IC, variable_setup.VAR, 0.0);
	V	= vector_2D(N_IC, variable_setup.VAR, 0.0);
	vector<CFD_mixture_data>	variable_data(N_IC);


	int ne	= variable_setup.NS+variable_setup.ND;

#pragma omp parallel for
	for (int n	= 0; n <= N_IC-1; n++)
	{
		for (int s = 0; s <= variable_setup.NS-1; s++)	V[n][s]						= rho_s[n][s];
		for (int k = 0; k <= variable_setup.ND-1; k++)	V[n][variable_setup.NS+k]	= u[n][k];

		V[n][ne]	= T[n];
		switch (variable_setup.energy_flag)
		{
		case 0:
			variable_data[n].calculate_data(variable_setup.NS, 0, V[n], species.Cv_t, species.Cv_r, species.R, species.M);
			break;

		case 1:
			V[n][ne+1]	= Te[n];
			variable_data[n].calculate_data(variable_setup.NS, 0, V[n], species.Cv_t, species.Cv_r, species.R, species.M);
			break;

		case 2:
			V[n][ne+1]	= Tv[n];
			variable_data[n].calculate_data(variable_setup.NS, 0, V[n], species.Cv_t, species.Cv_r, species.R, species.M);
			break;

		case 3:
			V[n][ne+1]	= Tv[n];
			V[n][ne+2]	= Te[n];
			variable_data[n].calculate_data(variable_setup.NS, 0, V[n], species.Cv_t, species.Cv_r, species.R, species.M);
			break;

		case 4:
			V[n][ne+1]	= Tr[n];
			variable_data[n].calculate_data(variable_setup.NS, 1, V[n], species.Cv_t, species.Cv_r, species.R, species.M);
			break;

		case 5:
			V[n][ne+1]	= Tr[n];
			V[n][ne+2]	= Te[n];
			variable_data[n].calculate_data(variable_setup.NS, 1, V[n], species.Cv_t, species.Cv_r, species.R, species.M);
			break;

		case 6:
			V[n][ne+1]	= Tr[n];
			V[n][ne+2]	= Tv[n];
			variable_data[n].calculate_data(variable_setup.NS, 1, V[n], species.Cv_t, species.Cv_r, species.R, species.M);
			break;

		case 7:
			V[n][ne+1]	= Tr[n];
			V[n][ne+2]	= Tv[n];
			V[n][ne+3]	= Te[n];
			variable_data[n].calculate_data(variable_setup.NS, 1, V[n], species.Cv_t, species.Cv_r, species.R, species.M);
			break;
		}
	}

	for (int n	= 0; n <= N_IC-1; n++)	CFD_V_to_Q(V[n], Q[n], variable_setup, variable_data[n], species);
}



void Calculate_IC_multi_fluid(vector<CFD_variable_setup_ver2>	&variable_setup,
								vector< vector < vector < double > > > &rho_s, vector< vector < vector < double > > > &u,
								vector< vector <double> > &T, vector< vector <double> > &Tr, vector< vector <double> > Tv, vector< vector <double> > Te,
								vector< vector< vector<double> > > &Q, vector< vector< vector<double> > > &V, vector<SPECIES_DATA_BASIC>	&species,
								int NUM_FLUID)
{
#pragma omp parallel for num_threads(NUM_FLUID)
	for (int f = 0; f <= NUM_FLUID-1; f++)
	{
		Calculate_IC(variable_setup[f],	rho_s[f], u[f], T[f], Tr[f], Tv[f], Te[f], Q[f], V[f], species[f]);
	}
}


void Calculate_IC_multi_fluid_ver2(vector<SOL_CFD> &Solution,
									vector< vector < vector < double > > > &rho_s, vector< vector < vector < double > > > &u,
									vector< vector <double> > &T, vector< vector <double> > &Tr, vector< vector <double> > &Tv, vector< vector <double> > &Te,
									vector< vector< vector<double> > > &Q, vector< vector< vector<double> > > &V, vector<SPECIES_DATA_BASIC>	&species,
									int NUM_FLUID)
{
#pragma omp parallel for num_threads(NUM_FLUID)
	for (int f = 0; f <= NUM_FLUID-1; f++)
	{
		Calculate_IC(Solution[f].setup,	rho_s[f], u[f], T[f], Tr[f], Tv[f], Te[f], Q[f], V[f], species[f]);
	}
}





void Assign_IC (CFD_variable_setup_ver2	variable_setup, vector< vector<double> > &Q, vector< vector<double> > &V,
				SOL_CLASS_BASIC	&Cell_data_Q, SOL_CLASS_BASIC	&Cell_data_V, SPECIES_DATA_BASIC &species, int NCM, int ini_method, int NT)
{

#pragma omp parallel for num_threads(NT)
	for (int c = 0; c <= NCM-1; c++)
	{
		for (int s = 0; s <= variable_setup.VAR-1; s++)
		{
			Cell_data_Q.data_ptr[c]->data[s]	= Q[ini_method][s];
			Cell_data_V.data_ptr[c]->data[s]	= V[ini_method][s];
		}
	}
}


void Assign_IC_ver2 (SOL_CFD &Solution, vector< vector<double> > &Q, vector< vector<double> > &V, int NCM, int ini_method, int NT)
{

#pragma omp parallel for num_threads(NT)
	for (int c = 0; c <= NCM-1; c++)
	{
		for (int s = 0; s <= Solution.setup.VAR; s++)
		{
			Solution.Qc[c][s]	= Q[ini_method][s];
			Solution.Vc[c][s]	= V[ini_method][s];
		}
	}
}

void Initialize_flows_multi_fluids(vector<SOL_CFD> &Solution,
									vector< vector< vector<double> > > &Q,
									vector< vector< vector<double> > > &V,
									vector<SPECIES_DATA_BASIC>	&species,
									int NCM, int ini_method, int NUM_FLUID, int NT)
{

	string line;
	ifstream restart_file;
	restart_file.open("restart.dat");

	int iter;
	int	nfluid;
	int type;
	int ncm;

	if (restart_file.is_open())
	{
		getline(restart_file, line);	data_read::read_line(line);
		data_read::get_data_from_string<int>(line,		"[ITERATION]: ", 0, iter);

		getline(restart_file, line);	data_read::read_line(line);
		data_read::get_data_from_string<int>(line,		"[NUMBER OF FLUID]:", 0, nfluid);

		getline(restart_file, line);	data_read::read_line(line);
		data_read::get_data_from_string<int>(line,		"[VARIABLES TYPE]:", 0, type);

		getline(restart_file, line);	data_read::read_line(line);
		data_read::get_data_from_string<int>(line,		"[NCM]:", 0, ncm);

		if (ncm != NCM)
		{
			cout << "PROBLEM IN RESTART DATA. TOTAL NUMBER OF CELL DATA DOES NOT MATHCH WITH MESH INFOMATION DATA" << endl;
			program_error_type1("Please check restart file");
		}


		for (int f = 0; f <= NUM_FLUID-1; f++)
		{

			switch (Solution[f].setup.energy_flag)
			{
			case 0:
				for (int c = 0; c <= NCM-1; c++)
				{
					for (int i = 0; i <= Solution[f].setup.VAR-1; i++)	restart_file >> Solution[f].Vc[c][i];
					Solution[f].mixture_data_c[c].calculate_data(Solution[f].setup.NS, 0, Solution[f].Vc[c], species[f].Cv_t, species[f].Cv_r, species[f].R, species[f].M);
				}
				break;

			case 1:
				for (int c = 0; c <= NCM-1; c++)
				{
					for (int i = 0; i <= Solution[f].setup.VAR-1; i++)	restart_file >> Solution[f].Vc[c][i];
					Solution[f].mixture_data_c[c].calculate_data(Solution[f].setup.NS, 0, Solution[f].Vc[c], species[f].Cv_t, species[f].Cv_r, species[f].R, species[f].M);
				}
				break;

			case 2:
				for (int c = 0; c <= NCM-1; c++)
				{
					for (int i = 0; i <= Solution[f].setup.VAR-1; i++)	restart_file >> Solution[f].Vc[c][i];
					Solution[f].mixture_data_c[c].calculate_data(Solution[f].setup.NS, 0, Solution[f].Vc[c], species[f].Cv_t, species[f].Cv_r, species[f].R, species[f].M);
				}
				break;

			case 3:
				for (int c = 0; c <= NCM-1; c++)
				{
					for (int i = 0; i <= Solution[f].setup.VAR-1; i++)	restart_file >> Solution[f].Vc[c][i];
					Solution[f].mixture_data_c[c].calculate_data(Solution[f].setup.NS, 0, Solution[f].Vc[c], species[f].Cv_t, species[f].Cv_r, species[f].R, species[f].M);
				}
				break;

			case 4:
				for (int c = 0; c <= NCM-1; c++)
				{
					for (int i = 0; i <= Solution[f].setup.VAR-1; i++)	restart_file >> Solution[f].Vc[c][i];
					Solution[f].mixture_data_c[c].calculate_data(Solution[f].setup.NS, 1, Solution[f].Vc[c], species[f].Cv_t, species[f].Cv_r, species[f].R, species[f].M);
				}
				break;

			case 5:
				for (int c = 0; c <= NCM-1; c++)
				{
					for (int i = 0; i <= Solution[f].setup.VAR-1; i++)	restart_file >> Solution[f].Vc[c][i];
					Solution[f].mixture_data_c[c].calculate_data(Solution[f].setup.NS, 1, Solution[f].Vc[c], species[f].Cv_t, species[f].Cv_r, species[f].R, species[f].M);
				}
				break;

			case 6:
				for (int c = 0; c <= NCM-1; c++)
				{
					for (int i = 0; i <= Solution[f].setup.VAR-1; i++)	restart_file >> Solution[f].Vc[c][i];
					Solution[f].mixture_data_c[c].calculate_data(Solution[f].setup.NS, 1, Solution[f].Vc[c], species[f].Cv_t, species[f].Cv_r, species[f].R, species[f].M);
				}
				break;

			case 7:
				for (int c = 0; c <= NCM-1; c++)
				{
					for (int i = 0; i <= Solution[f].setup.VAR-1; i++)	restart_file >> Solution[f].Vc[c][i];
					Solution[f].mixture_data_c[c].calculate_data(Solution[f].setup.NS, 1, Solution[f].Vc[c], species[f].Cv_t, species[f].Cv_r, species[f].R, species[f].M);
				}
				break;
			}

			for (int c = 0; c <= NCM-1; c++)	CFD_V_to_Q(Solution[f].Vc[c], Solution[f].Qc[c], Solution[f].setup, Solution[f].mixture_data_c[c], species[f]);
		}
	}
	else
	{
#pragma omp parallel for num_threads(NUM_FLUID)
		for (int f = 0; f <= NUM_FLUID-1; f++)	Assign_IC_ver2 (Solution[f], Q[f], V[f], NCM, ini_method, NT);
	}
}

