/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 20, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Jacobian_source_chemistry.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_CFD_Variable_change.hpp"
#include "../include/OP2A_CFD_reconstruct.hpp"
#include "../include/OP2A_CFD_Flux.hpp"
#include "../include/OP2A_CFD_Jacobian.hpp"


#include "../../UTIL/include/OP2A_utilities.hpp"
#include "Problem/include/OP2A_setup_problem.hpp"
#include "../../MATRIX/include/OPPA_Matrix.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../GRID/include/OP2A_grid.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"

//	Jacobian of Forward reactions
void CFD_Calculate_dRfb_dQ(double kfb, double Rfb, vector<double> &dkfb_dQ, vector<double> &alpha_j, vector<double> &rho_j, int NS, int VAR, vector<double> &dRfb_dQ)
{
	if (kfb != 0.0)
	{
		for (int i = 0;		i <= NS-1; 	i++)	if (alpha_j[i] != 0.0)	dRfb_dQ[i]	= dkfb_dQ[i]/kfb + rho_j[i]/alpha_j[i];
		for (int i = NS; 	i <= VAR-1; i++)	dRfb_dQ[i]	= dkfb_dQ[i]/kfb;

		for (int i = 0; 	i <= VAR-1; i++)	dRfb_dQ[i]	= Rfb * dRfb_dQ[i];
	}
	else
	{
		for (int i = 0; 	i <= VAR-1; i++)	dRfb_dQ[i]	= 0.0;
	}
}

void CFD_Calculate_dRfb_dQ_2(double kfb, double Rfb, vector<double> &dkfb_dQ, vector<double> &alpha_j, vector<double> &rho_j, int NS1, int ND1, int NE1, int VAR, vector<double> &dRfb_dQ)
{
	int start_second	= NS1 + ND1 + NE1;

	for (int i = 0;	i <= VAR-1; i++)	dRfb_dQ[i]	= dkfb_dQ[i]/kfb;

	for (int i = 0;	i <= NS1-1; i++)	dRfb_dQ[i]	+= rho_j[i]/alpha_j[i];
	dRfb_dQ[start_second]	+= rho_j[NS1]/alpha_j[NS1];

	for (int i = 0; 	i <= VAR-1; i++)	dRfb_dQ[i]	= Rfb * dRfb_dQ[i];
}


// Jacobian of forward reaction rate
void CFD_Calculate_dkf_dQ(double kf, double Tc, double nu_k, double theta_k, vector<double> &dTc_dQ, int NVAR, vector<double>	&dkf_dQ)
{
	double temp	= kf * (nu_k/Tc + theta_k/Tc/Tc);

	for (int i = 0; i <= NVAR-1; i++)	dkf_dQ[i]	= temp * dTc_dQ[i];
}


void CFD_Calculate_dkb_dQ(double kb, double Tbc, double nu_k, double theta_k, double dKeq_dT_over_Keq, vector<double> &dTbc_dQ, int NVAR, vector<double>	&dkb_dQ)
{
	double temp;

	temp	= kb * ((nu_k/Tbc + theta_k/Tbc/Tbc) + dKeq_dT_over_Keq);
	for (int i = 0; i <= NVAR-1; i++)	dkb_dQ[i]	= temp * dTbc_dQ[i];
}







void CFD_Jacobian_chemistry_single_fluid(GRID_CLASS &grid, SOL_CFD &Solution, vector< vector<double> > &kf, vector< vector<double> > &Rf, vector< vector<double> > &kb, vector< vector<double> > &Rb, REACTION_DATA_ver2 &reactions, vector< vector< vector<double> > > &dS_chem_dQ, bool is_axis, int nt)
{
	int ne	= Solution.setup.NS	+ Solution.setup.ND;

	double flag_axis;
	if (is_axis	== true)	flag_axis	= 1.0;
	else					flag_axis	= 0.0;

#pragma omp parallel for num_threads(nt)
	for (int c = 0; c <= grid.NCM-1; c++)
	{
		double Volume	= grid.cells.data_ptr[c]->S * (1.0 + flag_axis*(grid.cells.data_ptr[c]->x[1]- 1.0));

		double T_tra	= Solution.Vc[c][ne+Solution.setup.ID_T[TRA]];
		double T_rot	= Solution.Vc[c][ne+Solution.setup.ID_T[ROT]];
		double T_vib	= Solution.Vc[c][ne+Solution.setup.ID_T[VIB]];
		double T_e		= Solution.Vc[c][ne+Solution.setup.ID_T[ELE]];

		double n_mix = 0.0;
		for (int s = 0; s <= Solution.setup.NS-1; s++)	n_mix	+= Solution.Qc[c][s] / reactions.species_data[s].basic_data.m;

		// Initialize Jacobian of chemistry source terms
		for (int s = 0; s <= Solution.setup.NS-1; s++)
			for (int i = 0; i <= Solution.setup.VAR-1; i++) dS_chem_dQ[c][s][i]	= 0.0;


		for (int k = 0; k <= reactions.NR-1; k++)
		{
			double Tcf	= pow(T_tra, reactions.reaction_k[k].temperature_coeff_f[TRA])	* pow(T_rot, reactions.reaction_k[k].temperature_coeff_f[ROT]) * pow(T_vib, reactions.reaction_k[k].temperature_coeff_f[VIB]) * pow(T_e, reactions.reaction_k[k].temperature_coeff_f[ELE]);
			double Tcb	= pow(T_tra, reactions.reaction_k[k].temperature_coeff_b[TRA])	* pow(T_rot, reactions.reaction_k[k].temperature_coeff_b[ROT]) * pow(T_vib, reactions.reaction_k[k].temperature_coeff_b[VIB]) * pow(T_e, reactions.reaction_k[k].temperature_coeff_b[ELE]);

			// 1. Calculate dTc_dQ
			vector <double>	dTcf_dQ(Solution.setup.VAR, 0.0);
			vector <double>	dTcb_dQ(Solution.setup.VAR, 0.0);

			for (int i = 0; i <= Solution.setup.VAR-1; i++)
			{
				dTcf_dQ[i]	= Tcf * (reactions.reaction_k[k].temperature_coeff_f[0]/T_tra*Solution.dT_dQ[c][Solution.setup.ID_T[TRA]][i]
				          	       + reactions.reaction_k[k].temperature_coeff_f[1]/T_rot*Solution.dT_dQ[c][Solution.setup.ID_T[ROT]][i]
								   + reactions.reaction_k[k].temperature_coeff_f[2]/T_rot*Solution.dT_dQ[c][Solution.setup.ID_T[VIB]][i]
								   + reactions.reaction_k[k].temperature_coeff_f[3]/T_rot*Solution.dT_dQ[c][Solution.setup.ID_T[ELE]][i]);

				dTcb_dQ[i]	= Tcb * (reactions.reaction_k[k].temperature_coeff_b[0]/T_tra*Solution.dT_dQ[c][Solution.setup.ID_T[TRA]][i]
								   + reactions.reaction_k[k].temperature_coeff_b[1]/T_rot*Solution.dT_dQ[c][Solution.setup.ID_T[ROT]][i]
								   + reactions.reaction_k[k].temperature_coeff_b[2]/T_rot*Solution.dT_dQ[c][Solution.setup.ID_T[VIB]][i]
								   + reactions.reaction_k[k].temperature_coeff_b[3]/T_rot*Solution.dT_dQ[c][Solution.setup.ID_T[ELE]][i]);

			}


			// 2. Calcualte dkf_dQ
			vector <double>	dkf_dQ(Solution.setup.VAR, 0.0);
			vector <double>	dkb_dQ(Solution.setup.VAR, 0.0);

			double dKeq_dT_over_Keq	= reactions.reaction_k[k].Keq_data.calculate_dKeq_dT_over_Keq(Tcb, n_mix);
			CFD_Calculate_dkf_dQ(kf[c][k], Tcf, reactions.reaction_k[k].kf_coeff[1], reactions.reaction_k[k].kf_coeff[2], dTcf_dQ, Solution.setup.VAR, dkf_dQ);
			CFD_Calculate_dkb_dQ(kb[c][k], Tcb, reactions.reaction_k[k].kf_coeff[1], reactions.reaction_k[k].kf_coeff[2], dKeq_dT_over_Keq, dTcb_dQ, Solution.setup.VAR, dkb_dQ);

			// 3. Calculate dRf_dQ and dRb_dQ
			vector <double>	dRf_dQ(Solution.setup.VAR, 0.0);
			vector <double>	dRb_dQ(Solution.setup.VAR, 0.0);

			CFD_Calculate_dRfb_dQ(kf[c][k], Rf[c][k], dkf_dQ, reactions.reaction_k[k].Reactant_coeff, Solution.Vc[c], Solution.setup.NS, Solution.setup.VAR, dRf_dQ);
			CFD_Calculate_dRfb_dQ(kb[c][k], Rb[c][k], dkb_dQ, reactions.reaction_k[k].Product_coeff,  Solution.Vc[c], Solution.setup.NS, Solution.setup.VAR, dRb_dQ);



			for (int i = 0; i <= Solution.setup.VAR-1; i++)
			{
				for (int s = 0; s <= Solution.setup.NS-1; s++)
				{
					dS_chem_dQ[c][s][i]	+= Volume * reactions.species_data[s].basic_data.M * (reactions.reaction_k[k].Product_coeff[s] - reactions.reaction_k[k].Reactant_coeff[s]) * (dRf_dQ[i] - dRb_dQ[i]);
				}
			}
		}
	}
}













void CFD_Jacobian_dT_dQ_update_multi_fluid(int NCM, vector<CFD_variable_setup_ver2> &setup, vector<SOL_CFD> &Solutions, vector<SPECIES_DATA_BASIC> &species, REACTION_DATA_ver2 &reactions, vector< vector< vector<double> > > &dT_dQ, int num_fluid)
{
	int VAR = 0;
	for (int f = 0; f <= num_fluid-1; f++)	VAR += setup[f].VAR;

	dT_dQ.resize(NCM);
	for (int c = 0; c <= NCM-1; c++)	dT_dQ[c]	= vector_2D(4, VAR, 0.0);


	switch (num_fluid)
	{
	case 1:
		for (int c = 0; c <= NCM-1; c++)
		{
			for (int k = 0; k <= VAR-1; k++)
			{
				dT_dQ[c][TRA][k]	= Solutions[0].dT_dQ[c][0][k];
				dT_dQ[c][ROT][k]	= Solutions[0].dT_dQ[c][setup[0].ID_T[ROT]][k];
				dT_dQ[c][VIB][k]	= Solutions[0].dT_dQ[c][setup[0].ID_T[VIB]][k];
				dT_dQ[c][ELE][k]	= Solutions[0].dT_dQ[c][setup[0].ID_T[ELE]][k];
			}
		}
		break;

	case 2:
		for (int c = 0; c <= NCM-1; c++)
		{
			for (int k = 0; k <= setup[0].VAR-1; k++)
			{
				dT_dQ[c][TRA][k]	= Solutions[0].dT_dQ[c][0][k];
				dT_dQ[c][ROT][k]	= Solutions[0].dT_dQ[c][setup[0].ID_T[ROT]][k];
				dT_dQ[c][VIB][k]	= Solutions[0].dT_dQ[c][setup[0].ID_T[VIB]][k];
			}

			for (int k = 0; k <= setup[1].VAR-1; k++)
			{
				dT_dQ[c][ELE][setup[0].VAR + k]	= Solutions[1].dT_dQ[c][0][k];
			}
		}
		break;
	}
}

vector<double> CFD_Calculate_dTc_dQ(double a, double b, double c, double d, double T_tra, double T_rot, double T_vib, double T_e, vector<double> dT_tra_dQ, vector<double> dT_rot_dQ, vector<double> dT_vib_dQ, vector<double> dT_e_dQ, int NVAR)
{
	vector<double>	dTc_dQ(NVAR, 0.0);
	double	Tc	= pow(T_tra, a)	* pow(T_rot, b) * pow(T_vib, c) *pow(T_e, d);


	double a_T_tra	= a / T_tra;	if (T_tra == 0)	 a_T_tra 	= 0.0;
	double b_T_rot	= b / T_rot;	if (T_rot == 0)	 b_T_rot 	= 0.0;
	double c_T_vib	= c / T_vib;	if (T_vib == 0)	 c_T_vib 	= 0.0;
	double d_T_e	= d / T_e;		if (T_e == 0)	 d_T_e 		= 0.0;


#pragma omp parallel for
	for (int i = 0; i <= NVAR-1; i++)
	{
		dTc_dQ[i]	= Tc * (a_T_tra*dT_tra_dQ[i] + b_T_rot*dT_rot_dQ[i] + c_T_vib*dT_vib_dQ[i] + d_T_e*dT_e_dQ[i]);
	}

	return (dTc_dQ);
}














/*

void CFD_Jacobian_chemistry(CFD_variable_setup_ver2 setup, vector <double *> &rho_s_ALL, vector<vector <double *> > &Ts_ALL,
							vector< vector<double> > &dT_dQ,
							vector<double> &kf, vector<double> &Rf, vector<double> &kb, vector<double> &Rb,
							SPECIES_DATA_BASIC	&species, REACTION_DATA_ver2 &reactions, vector< vector<double> > &dS_chem_dQ)
{
	vector< vector<double> > dkf_dQ		= vector_2D(reactions.NR, setup.VAR, 0.0);
	vector< vector<double> > dRfkf_dQ	= vector_2D(reactions.NR, setup.VAR, 0.0);
	vector< vector<double> > dRf_dQ		= vector_2D(reactions.NR, setup.VAR, 0.0);
	vector< vector<double> > dRb_dQ		= vector_2D(reactions.NR, setup.VAR,  0.0);



	for (int k = 0; k <= reactions.NR-1; k++)
	{
		double T_tra	= 0.0;
		double T_rot	= 0.0;
		double T_vib	= 0.0;
		double T_e		= 0.0;

		for (int s = 0; s <= setup.NS-1; s++)
		{
			T_tra	+= *Ts_ALL[s][TRA] * reactions.Xf_all[k][s];
			T_rot	+= *Ts_ALL[s][ROT] * reactions.Xf_molecule[k][s];
			T_vib	+= *Ts_ALL[s][VIB] * reactions.Xf_molecule[k][s];
			T_e		+= *Ts_ALL[s][ELE] * reactions.Xf_electron[k][s];
		}

		if (T_rot 	== 0.0)	T_rot 	= 1.0;
		if (T_vib 	== 0.0)	T_vib 	= 1.0;
		if (T_e 	== 0.0)	T_e 	= 1.0;

		double Tcf	= pow(T_tra, reactions.reaction_k[k].temperature_coeff_f[TRA])	* pow(T_rot, reactions.reaction_k[k].temperature_coeff_f[ROT]) * pow(T_vib, reactions.reaction_k[k].temperature_coeff_f[VIB]) * pow(T_e, reactions.reaction_k[k].temperature_coeff_f[ELE]);
		double Tcb	= pow(T_tra, reactions.reaction_k[k].temperature_coeff_b[TRA])	* pow(T_rot, reactions.reaction_k[k].temperature_coeff_b[ROT]) * pow(T_vib, reactions.reaction_k[k].temperature_coeff_b[VIB]) * pow(T_e, reactions.reaction_k[k].temperature_coeff_b[ELE]);

		vector <double>	dTcf_dQ	= CFD_Calculate_dTc_dQ(reactions.reaction_k[k].temperature_coeff_f[0], reactions.reaction_k[k].temperature_coeff_f[1], reactions.reaction_k[k].temperature_coeff_f[2], reactions.reaction_k[k].temperature_coeff_f[3], T_tra, T_rot, T_vib, T_e, dT_dQ[setup.ID_T[TRA]], dT_dQ[setup.ID_T[ROT]], dT_dQ[setup.ID_T[VIB]], dT_dQ[setup.ID_T[ELE]], setup.VAR);
		vector <double>	dTcb_dQ	= CFD_Calculate_dTc_dQ(reactions.reaction_k[k].temperature_coeff_b[0], reactions.reaction_k[k].temperature_coeff_b[1], reactions.reaction_k[k].temperature_coeff_b[2], reactions.reaction_k[k].temperature_coeff_b[3], T_tra, T_rot, T_vib, T_e, dT_dQ[setup.ID_T[TRA]], dT_dQ[setup.ID_T[ROT]], dT_dQ[setup.ID_T[VIB]], dT_dQ[setup.ID_T[ELE]], setup.VAR);

		// 1. Forward reaction
		double temp1	= kf * (reactions.reaction_k[k].kf_coeff[1]/Tcf + reactions.reaction_k[k].kf_coeff[2]/Tcf/Tcf);
		for (int i = 0; i <= setup.VAR-1; i++)	dkf_dQ[i]	= temp1 * dTcf_dQ[i];

		for (int s = 0; s <= setup.NS-1; s++)
		{
			if (reactions.reaction_k[k].Reactant_coeff[s] != 0)	dRfkf_dQ[s] = (reactions.reaction_k[k].Reactant_coeff[s]/(*rho_s_ALL[s])) * Rf[k];
		}

		double temp2	= Rf[k] / kf[k];
		for (int i = 0; i <= setup.VAR-1; i++)	dRf_dQ[i]	= temp1*dkf_dQ[i]	+ dRfkf_dQ[i];


		// 2. Backward reaction
		double Rbkb	= Rb[k] / (kb[k] * 1000.0);
	}




	for (int s = 0; s <= setup.NS; s++)
	{
		for (int i = 0; i <= setup.VAR-1; i++)	dS_chem_dQ[s][i]	= 0.0;

		for (int k = 0; k <= reactions.NR-1; k++)
		{
			double temp	= reactions.reaction_k[k].Product_coeff[s] - reactions.reaction_k[k].Reactant_coeff[s];
			for (int i = 0; i <= setup.VAR-1; i++)	dS_chem_dQ[s][i]	+= temp * (dRf_dQ[k][i] - dRb_dQ[k][i]);
		}

		for (int i = 0; i <= setup.VAR-1; i++)	dS_chem_dQ[s][i]	*= species.M[s];
	}
}

*/

