/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 1, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Flux_inviscid_AUSM.cpp
 * 			-  
 *  
 */



#include "../include/OP2A_CFD_Variable_change.hpp"
#include "../include/OP2A_CFD_Flux.hpp"

#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../../SETUP/include/OP2A_setup_problem.hpp"
#include "../../MATRIX/include/OPPA_Matrix.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../GRID/include/OP2A_grid.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"

void CFD_flux_1_2_AUSM_single(CFD_variable_setup_ver2 setup,
							vector<double>	&Q_L, vector<double> &U_L,  vector<double> &W_L, CFD_mixture_data &variable_data_L,
							vector<double>	&Q_R, vector<double> &U_R,  vector<double> &W_R, CFD_mixture_data &variable_data_R,
							vector<double>	&F_1_2, vector< vector<double> > &n, double S, SPECIES_DATA_BASIC	&species)
{
	double alpha = 3.0/16.0;
	double beta	 = 1.0/8.0;
	double Ku	 = 0.75;
	double Kp	 = 0.25;

	double M_inf = 10.0;
	double sigma = 1.0;

	int ne	= setup.NS + setup.ND;

	// Calculate Normal/Tangential/Tangential2 velocities
	vector<double>	Up_L(setup.ND, 0.0);
	vector<double>	Up_R(setup.ND, 0.0);

	for (int k1 = 0; k1 <= setup.ND-1; k1++)
	{
		for (int k = 0; k <= setup.ND-1; k++)
		{
			Up_L[k1]	+= U_L[setup.NS+k] * n[k1][k];
			Up_R[k1]	+= U_R[setup.NS+k] * n[k1][k];
		}
	}

	// Absolute values of velocities
	double Vabs_L, Vabs_R;
	Vabs_L	= 0.0;
	Vabs_R	= 0.0;
	for (int k = 0; k <= setup.ND-1; k++)
	{
		Vabs_L	+= pow(U_L[setup.NS+k], 2.0);
		Vabs_R	+= pow(U_R[setup.NS+k], 2.0);
	}

	// Enthalphy
	double Ht_L, Ht_R;
	Ht_L	= (Q_L[ne] + W_L[ne])/variable_data_L.rho;//	- KE_L;
	Ht_R	= (Q_R[ne] + W_R[ne])/variable_data_R.rho;//	- KE_R;





	// 1. Calculate a_L and a_R
	vector < vector <double> >	dTp_dQ = vector_2D(setup.NE, setup.VAR, 0.0);
	vector < vector <double> >	dTm_dQ = vector_2D(setup.NE, setup.VAR, 0.0);
	CFD_Calculate_dT_dQ(setup, Q_L, U_L, variable_data_L, dTp_dQ, species);
	CFD_Calculate_dT_dQ(setup, Q_R, U_R, variable_data_R, dTm_dQ, species);

	vector<double>	dpp_dQ(setup.VAR, 0.0);
	vector<double>	dpm_dQ(setup.VAR, 0.0);
	CFD_Calculate_dp_dQ(setup, U_L, variable_data_L, dTp_dQ, dpp_dQ,	species);
	CFD_Calculate_dp_dQ(setup, U_R, variable_data_R, dTm_dQ, dpm_dQ,	species);

	double a2_L		= CFD_Calculate_a2(setup, Q_L, W_L, variable_data_L, dpp_dQ, species);
	double a2_R		= CFD_Calculate_a2(setup, Q_R, W_R, variable_data_R, dpm_dQ, species);

	double a_L;
	if (a2_L >= 0.0)	a_L	= sqrt(a2_L);
	else			program_error_type1("Error negative a2 LEFT [AUSM SCHEME]: Please check dp and dT calculation!!");

	double a_R;
	if (a2_R >= 0.0)	a_R	= sqrt(a2_R);
	else			program_error_type1("Error negative a2 RIGHT [AUSM SCHEME]: Please check dp and dT calculation!!");

	double c_bar	= (a_L + a_R) / 2.0;


	double Mp, Mm;
	Mp	= Up_L[0] / c_bar;
	Mm	= Up_R[0] / c_bar;

	double g 		= -fmax(fmin(Mp, 0.0), -1.0) * fmin(fmax(Mm, 0.0), 1.0);
	double Vn_mean	= (variable_data_L.rho*fabs(Up_L[0]) + variable_data_R.rho*fabs(Up_R[0])) / (variable_data_L.rho + variable_data_R.rho);

	double M_hat;
	double temp1;
	temp1 	= sqrt((Vabs_L + Vabs_R)/2.0);
	M_hat	= fmin(1.0, temp1/c_bar);

	double Xai	= pow((1.0 - M_hat), 2.0);

	double delta_p		= W_L[ne]	- W_R[ne];

	double Vn_p		= (1.0 - g)*Vn_mean + g*fabs(Up_L[0]);
	double Vn_m		= (1.0 - g)*Vn_mean + g*fabs(Up_R[0]);

	double m_dot	= 0.5 * (variable_data_L.rho*(Up_L[0] + Vn_p) + variable_data_R.rho*(Up_R[0] - Vn_m) - Xai/c_bar*delta_p);


	double beta_L;
	if (fabs(Mp) >= 1.0)	beta_L	= 0.5 * (1.0 + fsgn(Mp));
	else					beta_L	= 0.25 * (2.0 - Mp) * pow((Mp+1.0), 2.0);

	double beta_R;
	if (fabs(Mm) >= 1.0)	beta_R	= 0.5 * (1.0 - fsgn(Mm));
	else					beta_R	= 0.25 * (2.0 + Mm) * pow((Mm-1.0), 2.0);


	double c_star2_L	= 2.0*(variable_data_L.gamma - 1.0)/(variable_data_L.gamma + 1.0) * Ht_L;
	double c_star2_R	= 2.0*(variable_data_R.gamma - 1.0)/(variable_data_R.gamma + 1.0) * Ht_R;
	double c_L			= c_star2_L / fmax(sqrt(c_star2_L), fabs(Up_L[0]));
	double c_R			= c_star2_R / fmax(sqrt(c_star2_R), fabs(Up_R[0]));
	double c_1_2		= fmin(c_L, c_R);


	double p_1_2;
	double rho_bar 	= (variable_data_L.rho + variable_data_R.rho) / 2.0;
	p_1_2	= (W_L[ne]	+ W_R[ne])/2.0 + (beta_L - beta_R)/2.0*(W_L[ne] - W_R[ne]);// + sqrt(0.5*(Vabs_L + Vabs_R))*(beta_L+beta_R-1.0)*rho_bar*c_1_2;


	// 7 Phi_L and Phi_R
	vector<double>	Phi_L(setup.VAR, 0.0);
	vector<double>	Phi_R(setup.VAR, 0.0);

	for (int i = 0; i <= setup.VAR-1; i++)
	{
		Phi_L[i]	= Q_L[i] / variable_data_L.rho;
		Phi_R[i]	= Q_R[i] / variable_data_R.rho;
	}

	Phi_L[ne]	+= W_L[ne] / variable_data_L.rho;
	Phi_R[ne]	+= W_R[ne] / variable_data_R.rho;

	if (fabs(p_1_2) > 1.0)
	{
		int aa = 0;
	}
	for (int i = 0; i <= setup.VAR-1; i++)	F_1_2[i]			= (m_dot + fabs(m_dot))/2.0*Phi_L[i] + (m_dot - fabs(m_dot))/2.0*Phi_R[i];
	for (int k = 0; k <= setup.ND-1; k++)	F_1_2[setup.NS+k]	+= p_1_2 * n[0][k];
}
