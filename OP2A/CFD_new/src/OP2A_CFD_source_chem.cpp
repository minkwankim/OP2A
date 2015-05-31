/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 16, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_source_chem.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_CFD_Variable_change.hpp"
#include "../include/OP2A_CFD_reconstruct.hpp"
#include "../include/OP2A_CFD_Flux.hpp"


#include "../../UTIL/include/OP2A_utilities.hpp"
#include "Problem/include/OP2A_setup_problem.hpp"
#include "../../MATRIX/include/OPPA_Matrix.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../GRID/include/OP2A_grid.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"





void CFD_source_chemical_NONEQ_ver2(vector <double *> &rho_s_ALL, vector<vector <double *> > &Ts_ALL, REACTION_DATA_ver2 &reactions, vector<double> &kf, vector<double> &Rf, vector<double> &kb, vector<double> &Rb, vector<double *> &S_chem, int nt)
{
	int NR	= reactions.NR;
	int NS	= reactions.NS;

	vector<double> Tcf(NR, 0.0);
	vector<double> Tcb(NR, 0.0);


	// 1. Calculate forward and backward control temperatures
	for(int k = 0; k <= NR-1; k++)
	{
		// Forward control temperature
		double T_tra	= 0.0;
		double T_rot	= 0.0;
		double T_vib	= 0.0;
		double T_e		= 0.0;

		for (int s = 0; s <= NS-1; s++)
		{
			T_tra	+= *Ts_ALL[s][TRA] * reactions.Xf_all[k][s];
			T_rot	+= *Ts_ALL[s][ROT] * reactions.Xf_molecule[k][s];
			T_vib	+= *Ts_ALL[s][VIB] * reactions.Xf_molecule[k][s];
			T_e		+= *Ts_ALL[s][ELE] * reactions.Xf_electron[k][s];
		}

		if (T_rot 	== 0.0)	T_rot 	= 1.0;
		if (T_vib 	== 0.0)	T_vib 	= 1.0;
		if (T_e 	== 0.0)	T_e 	= 1.0;

		Tcf[k]	= pow(T_tra, reactions.reaction_k[k].temperature_coeff_f[TRA])	* pow(T_rot, reactions.reaction_k[k].temperature_coeff_f[ROT]) * pow(T_vib, reactions.reaction_k[k].temperature_coeff_f[VIB]) * pow(T_e, reactions.reaction_k[k].temperature_coeff_f[ELE]);

		// BAckward control temperature
		T_tra	= 0.0;
		T_rot	= 0.0;
		T_vib	= 0.0;
		T_e		= 0.0;

		for (int s = 0; s <= NS-1; s++)
		{
			T_tra	+= *Ts_ALL[s][TRA] * reactions.Xb_all[k][s];
			T_rot	+= *Ts_ALL[s][ROT] * reactions.Xb_molecule[k][s];
			T_vib	+= *Ts_ALL[s][VIB] * reactions.Xb_molecule[k][s];
			T_e		+= *Ts_ALL[s][ELE] * reactions.Xb_electron[k][s];
		}

		if (T_rot 	== 0.0)	T_rot 	= 1.0;
		if (T_vib 	== 0.0)	T_vib 	= 1.0;
		if (T_e 	== 0.0)	T_e 	= 1.0;

		Tcb[k]	= pow(T_tra, reactions.reaction_k[k].temperature_coeff_b[TRA])	* pow(T_rot, reactions.reaction_k[k].temperature_coeff_b[ROT]) * pow(T_vib, reactions.reaction_k[k].temperature_coeff_b[VIB]) * pow(T_e, reactions.reaction_k[k].temperature_coeff_b[ELE]);
	}


	// 2. Calculate forward reaction rates
	for(int k = 0; k <= NR-1; k++)
	{
		 kf[k]	= cal_kf(reactions.reaction_k[k], Tcf[k]);
	}



	// 3. Calculate backward reaction rates
	double n_mix = 0.0;
	for (int s = 0; s <= reactions.NS-1; s++)
		n_mix	+= (*rho_s_ALL[s]) / reactions.species_data[s].basic_data.m;

	for(int k = 0; k <= NR-1; k++)
	{
		 kb[k]	= cal_kb(reactions.reaction_k[k], n_mix, Tcb[k]);
	}


	// 4. Reaction rates
	vector<double>	rhos_M(NS, 0.0);
	for(int s = 0; s <= NS-1; s++)	rhos_M[s]	= *rho_s_ALL[s] / reactions.species_data[s].basic_data.M * 1.0e-3;	// [mol/cm^3]


	for(int k = 0; k <= NR-1; k++)
	{
		Rf[k]	= 1000.0 * kf[k];
		Rb[k]	= 1000.0 * kb[k];

		for (int s = 0; s <= NS-1; s++)
		{
			Rf[k]	*= pow(rhos_M[s], reactions.reaction_k[k].Reactant_coeff[s]);
			Rb[k]	*= pow(rhos_M[s], reactions.reaction_k[k].Product_coeff[s]);
		}
	}

	vector<double>	Rf_m_Rb(NR, 0.0);
	for(int k = 0; k <= NR-1; k++)	Rf_m_Rb[k]	= Rf[k] - Rb[k];


	// 5. Calculate Chemistry source term
	for (int s = 0; s <= NS-1; s++)
	{
		*S_chem[s]	= 0.0;

		for (int k = 0; k <= NR-1; k++)
		{
			double beta_m_alpha	= reactions.reaction_k[k].Product_coeff[s]	- reactions.reaction_k[k].Reactant_coeff[s];
			*S_chem[s]	+= beta_m_alpha * Rf_m_Rb[k];
		}

		*S_chem[s]	= reactions.species_data[s].basic_data.M * (*S_chem[s]);
	}
}


void CFD_source_chemical_NONEQ_ver1(vector <double *> &rho_s_ALL, vector<vector <double *> > &Ts_ALL, REACTION_DATA_ver2 &reactions, vector<double *> &S_chem, int nt)
{
	int NR	= reactions.NR;
	int NS	= reactions.NS;

	vector<double> Tcf(NR, 0.0);
	vector<double> Tcb(NR, 0.0);

	vector <double> Rf(NR, 0.0);
	vector <double> Rb(NR, 0.0);


	// 1. Calculate forward and backward control temperatures
//#pragma omp parallel for num_threads(nt)
	for(int k = 0; k <= NR-1; k++)
	{
		// Forward control temperature
		double T_tra	= 0.0;
		double T_rot	= 0.0;
		double T_vib	= 0.0;
		double T_e		= 0.0;

		for (int s = 0; s <= NS-1; s++)
		{
			T_tra	+= *Ts_ALL[s][TRA] * reactions.Xf_all[k][s];
			T_rot	+= *Ts_ALL[s][ROT] * reactions.Xf_molecule[k][s];
			T_vib	+= *Ts_ALL[s][VIB] * reactions.Xf_molecule[k][s];
			T_e		+= *Ts_ALL[s][ELE] * reactions.Xf_electron[k][s];
		}

		if (T_rot 	== 0.0)	T_rot 	= 1.0;
		if (T_vib 	== 0.0)	T_vib 	= 1.0;
		if (T_e 	== 0.0)	T_e 	= 1.0;

		Tcf[k]	= pow(T_tra, reactions.reaction_k[k].temperature_coeff_f[TRA])	* pow(T_rot, reactions.reaction_k[k].temperature_coeff_f[ROT]) * pow(T_vib, reactions.reaction_k[k].temperature_coeff_f[VIB]) * pow(T_e, reactions.reaction_k[k].temperature_coeff_f[ELE]);

		// BAckward control temperature
		T_tra	= 0.0;
		T_rot	= 0.0;
		T_vib	= 0.0;
		T_e		= 0.0;

		for (int s = 0; s <= NS-1; s++)
		{
			T_tra	+= *Ts_ALL[s][TRA] * reactions.Xb_all[k][s];
			T_rot	+= *Ts_ALL[s][ROT] * reactions.Xb_molecule[k][s];
			T_vib	+= *Ts_ALL[s][VIB] * reactions.Xb_molecule[k][s];
			T_e		+= *Ts_ALL[s][ELE] * reactions.Xb_electron[k][s];
		}

		if (T_rot 	== 0.0)	T_rot 	= 1.0;
		if (T_vib 	== 0.0)	T_vib 	= 1.0;
		if (T_e 	== 0.0)	T_e 	= 1.0;

		Tcb[k]	= pow(T_tra, reactions.reaction_k[k].temperature_coeff_b[TRA])	* pow(T_rot, reactions.reaction_k[k].temperature_coeff_b[ROT]) * pow(T_vib, reactions.reaction_k[k].temperature_coeff_b[VIB]) * pow(T_e, reactions.reaction_k[k].temperature_coeff_b[ELE]);
	}


	// 2. Calculate forward reaction rates
	vector<double>	kf(NR, 0.0);
//#pragma omp parallel for num_threads(nt)
	for(int k = 0; k <= NR-1; k++)
	{
		 kf[k]	= cal_kf(reactions.reaction_k[k], Tcf[k]);
	}



	// 3. Calculate backward reaction rates
	double n_mix = 0.0;
	for (int s = 0; s <= reactions.NS-1; s++)
		n_mix	+= (*rho_s_ALL[s]) / reactions.species_data[s].basic_data.m;

	vector<double>	kb(NR, 0.0);
//#pragma omp parallel for num_threads(nt)
	for(int k = 0; k <= NR-1; k++)
	{
		kb[k]	= cal_kb(reactions.reaction_k[k], n_mix, Tcb[k]);
	}


	// 4. Reaction rates
	vector<double>	rhos_M(NS, 0.0);
//#pragma omp parallel for num_threads(nt)
	for(int s = 0; s <= NS-1; s++)	rhos_M[s]	= *rho_s_ALL[s] / reactions.species_data[s].basic_data.M * 1.0e-3;	// [mol/cm^3]


//#pragma omp parallel for num_threads(nt)
	for(int k = 0; k <= NR-1; k++)
	{
		Rf[k]	= 1000.0 * kf[k];
		Rb[k]	= 1000.0 * kb[k];

		for (int s = 0; s <= NS-1; s++)
		{
			Rf[k]	*= pow(rhos_M[s], reactions.reaction_k[k].Reactant_coeff[s]);
			Rb[k]	*= pow(rhos_M[s], reactions.reaction_k[k].Product_coeff[s]);
		}
	}

	vector<double>	Rf_m_Rb(NR, 0.0);
//#pragma omp parallel for num_threads(nt)
	for(int k = 0; k <= NR-1; k++)	Rf_m_Rb[k]	= Rf[k] - Rb[k];


	// 5. Calculate Chemistry source term
	for (int s = 0; s <= NS-1; s++)
	{
		*S_chem[s]	= 0.0;

		for (int k = 0; k <= NR-1; k++)
		{
			double beta_m_alpha	= reactions.reaction_k[k].Product_coeff[s]	- reactions.reaction_k[k].Reactant_coeff[s];
			*S_chem[s]	+= beta_m_alpha * Rf_m_Rb[k];
		}

		*S_chem[s]	= reactions.species_data[s].basic_data.M * (*S_chem[s]);
	}
}








void CFD_source_chemical_NONEQ_ALL_explicit(int NCM, vector< vector <double *> > &rho_s_ALL, vector< vector<vector <double *> > > &Ts_ALL, REACTION_DATA_ver2 &reactions, vector< vector<double *> > &S_chem, int nt)
{
#pragma omp parallel for num_threads(nt)
	for (int c = 0; c <= NCM-1; c++)
	{
		CFD_source_chemical_NONEQ_ver1(rho_s_ALL[c], Ts_ALL[c], reactions, S_chem[c], nt);
	}
}

void CFD_source_chemical_NONEQ_ALL_implicit(int NCM, vector< vector <double *> > &rho_s_ALL, vector< vector<vector <double *> > > &Ts_ALL, REACTION_DATA_ver2 &reactions,
											vector<vector<double> > &kf, vector<vector<double> > &Rf, vector<vector<double> > &kb, vector<vector<double> > &Rb,
											vector< vector<double *> > &S_chem, int nt)
{
#pragma omp parallel for num_threads(nt)
	for (int c = 0; c <= NCM-1; c++)
	{
		CFD_source_chemical_NONEQ_ver2(rho_s_ALL[c], Ts_ALL[c], reactions, kf[c], Rf[c], kb[c], Rb[c], S_chem[c], nt);
	}
}



void CFD_residue_CHEM_NONEQ(GRID_CLASS &grid_data, SOL_CFD &Solutions, bool is_axis, int nt)
{
#pragma omp parallel for num_threads(nt)
	for (int c = 0; c <= grid_data.NCM-1; c++)
	{
		double S = 0.0;

		if (is_axis	== true)	S	= grid_data.cells.data_ptr[c]->S * fabs(grid_data.cells.data_ptr[c]->x[1]);
		else					S	= grid_data.cells.data_ptr[c]->S;

		for (int s = 0; s <= Solutions.setup.NS-1; s++)
		{
			double Source_term	= Solutions.S_source[c][s];
			Source_term	= Source_term * S;

			Solutions.Rn[c][s]	-= Source_term;
		}
	}
}
