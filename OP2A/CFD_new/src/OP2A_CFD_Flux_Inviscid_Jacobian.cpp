/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 15, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Flux_Inviscid_Jacobian.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_CFD_Variable_change.hpp"
#include "../include/OP2A_CFD_Flux.hpp"

#include "../../UTIL/include/OP2A_utilities.hpp"
#include "Problem/include/OP2A_setup_problem.hpp"
#include "../../MATRIX/include/OPPA_Matrix.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../GRID/include/OP2A_grid.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"



void CFD_Jacobian_inviscid_single(CFD_variable_setup_ver2 setup,	vector<double>	&Q, vector<double> &V, vector<double>	&W, CFD_mixture_data &variable_data,
									vector < vector<double> >	&dT_dQ, vector<double> &dp_dQ,
									vector < vector<double> >	&Jacob_inviscid,
									double signal, vector< vector<double> > &n, SPECIES_DATA_BASIC	&species,
									double dp, double dist_wall, double n_dot_wall, double x0, double eps0)
{
	// 1. Calculate Normal/Tangential/Tangential2 velocities
	vector<double>	Up(setup.ND, 0.0);

	for (int k1 = 0; k1 <= setup.ND-1; k1++)
	{
		for (int k = 0; k <= setup.ND-1; k++)
		{
			Up[k1]	+= V[setup.NS+k] * n[k1][k];
		}
	}


	// 2. Calculate d2p_dQ2
	vector < vector<double>	> 	d2p_dQ2 = vector_2D(setup.VAR, setup.VAR, 0.0);
	CFD_Calculate_d2p_dQ2(setup, V, variable_data, d2p_dQ2,	species);


	// 3. Calculate a
	double a2		= CFD_Calculate_a2(setup, Q, W, variable_data, dp_dQ, species);
	double two_a2	= 2.0 *a2;
	double a;
	if (a2 >= 0.0)	a	= sqrt(a2);
	else			program_error_type1("Error negative a2: Please check dp and dT calculation!!");


	// 4. Calculate da2_dQ and da_dQ
	vector<double> da2_dQ(setup.VAR, 0.0);
	vector<double> da_dQ(setup.VAR, 0.0);

	CFD_Calculate_da2_dQ(setup, Q, V, W, variable_data, dp_dQ, d2p_dQ2,	da2_dQ, species);
	CFD_Calculate_da_dQ(setup, Q, V, W, variable_data, da2_dQ, a, da_dQ, species);


	// 5. Calculate eigenvalues
	double lambda1, lambda2, lambda3;
	double eps;
	if (dist_wall < x0)	eps	= eps0 * n_dot_wall * (fabs(Up[0]) + a);
	else				eps	= eps0 * (fabs(Up[0]) + a);

	lambda1	= 0.5 * (Up[0] 		+ signal*sqrt(Up[0]*Up[0] 		+ eps*eps));
	lambda2	= 0.5 * ((Up[0]+a)	+ signal*sqrt(pow(Up[0]+a, 2.0) + eps*eps));
	lambda3	= 0.5 * ((Up[0]-a)	+ signal*sqrt(pow(Up[0]-a, 2.0) + eps*eps));

	vector <double> dU_dQ(setup.VAR, 0.0);
	double temp1 = Up[0]/variable_data.rho;

#pragma omp parallel for
	for (int s = 0; s <= setup.NS-1; s++)	dU_dQ[s]			= -temp1;

#pragma omp parallel for
	for (int k = 0; k <= setup.ND-1; k++)	dU_dQ[setup.NS+k]	= n[0][k]/variable_data.rho;

	vector<double>	dlambda1(setup.VAR, 0.0);
	vector<double>	dlambda2(setup.VAR, 0.0);
	vector<double>	dlambda3(setup.VAR, 0.0);

	temp1	= Up[0] / sqrt(Up[0]*Up[0] + eps*eps);
#pragma omp parallel for
	for (int i = 0; i <= setup.VAR-1; i++)	dlambda1[i]	= 0.5*(1.0 + signal*temp1)*dU_dQ[i];

	temp1	= (Up[0]+a) / sqrt((Up[0]+a)*(Up[0]+a) + eps*eps);
#pragma omp parallel for
	for (int i = 0; i <= setup.VAR-1; i++)	dlambda2[i]	= 0.5*(1.0 + signal*temp1) * (dU_dQ[i]+da_dQ[i]);

	temp1	= (Up[0]-a) / sqrt((Up[0]-a)*(Up[0]-a) + eps*eps);
#pragma omp parallel for
	for (int i = 0; i <= setup.VAR-1; i++)	dlambda3[i]	= 0.5*(1.0 + signal*temp1) * (dU_dQ[i]-da_dQ[i]);


	/*
	 * [PRE] Calculate predefined variables
	 */
	int ne = setup.NS + setup.ND;
	double A	= W[ne]/a2 * dp_dQ[ne];
	double B	= 0.5 * (variable_data.rho - A);
	vector <double> dA_dQ(setup.VAR, 0.0);
	vector <double> dB_dQ(setup.VAR, 0.0);
	vector <double> drho_dQ(setup.VAR, 0.0);

#pragma omp parallel for
	for (int i = 0; i <= setup.VAR-1;	i++)	dA_dQ[i]	= (dp_dQ[i]*a2 - da2_dQ[i]*W[ne])/a2/a2*dp_dQ[ne]	+ (W[ne]/a2)*d2p_dQ2[ne][i];

#pragma omp parallel for
	for (int i = 0; i <= setup.NS-1; 	i++)	drho_dQ[i]	= 1.0;

#pragma omp parallel for
	for (int i = 0; i <= setup.VAR-1;	i++)	dB_dQ[i]	= 0.5 * (drho_dQ[i] - dA_dQ[i]);

	vector <double> Qp(setup.VAR, 0.0);
#pragma omp parallel for
	for (int i = 0; i <= setup.VAR-1;	i++)	Qp[i]	= Q[i] / variable_data.rho;
	Qp[ne]	+= W[ne]/variable_data.rho;


	vector < vector <double> > dQp_dQ	= vector_2D(setup.VAR, setup.VAR, 0.0);
#pragma omp parallel for
	for (int r = 0; r <= setup.VAR-1;	r++)
	{
		for (int s = 0; s <= setup.VAR-1;	s++)
		{
			dQp_dQ[r][s]	= (Delta_fn(s,r)*variable_data.rho	- drho_dQ[s]*Q[r]) / pow(variable_data.rho, 2.0);
		}
	}

	for (int s = 0; s <= setup.VAR-1;	s++)	dQp_dQ[ne][s]	+= (dp_dQ[s]*variable_data.rho	- drho_dQ[s]*W[ne]) / pow(variable_data.rho, 2.0);

	vector <double> dUa_dQ(setup.VAR, 0.0);
#pragma omp parallel for
	for (int i = 0; i <= setup.VAR-1;	i++)	dUa_dQ[i]	= dU_dQ[i]*a	+ Up[0]*da_dQ[i];



	/*
	 * A. Species
	 */
#pragma omp parallel for
	for (int s = 0; s <= setup.NS-1; s++)
	{
		for (int i = 0; i <= setup.VAR-1; i++)
		{
			Jacob_inviscid[s][i] =  (A*dQp_dQ[s][i] + dA_dQ[i]*Qp[s])*lambda1	+ A*Qp[s]*dlambda1[i]
			                      + (B*dQp_dQ[s][i] + dB_dQ[i]*Qp[s])*(lambda2 + lambda3)	+ B*Qp[s]*(dlambda2[i] + dlambda3[i]);

		}
	}


	/*
	 * B. Momentum
	 */
#pragma omp parallel for
	for (int k = 0; k <= setup.ND-1; k++)
	{
		for (int i = 0; i <= setup.VAR-1; i++)
		{
			Jacob_inviscid[setup.NS+k][i] =  (A*dQp_dQ[setup.NS+k][i] + dA_dQ[i]*Qp[setup.NS+k])*lambda1	+ A*Qp[setup.NS+k]*dlambda1[i];
			Jacob_inviscid[setup.NS+k][i] += dB_dQ[i] * ((Qp[setup.NS+k] + a*n[0][k])*lambda2 + (Qp[setup.NS+k] - a*n[0][k])*lambda3);
			Jacob_inviscid[setup.NS+k][i] += B * ((dQp_dQ[setup.NS+k][i] + da_dQ[i]*n[0][k])*lambda2 + (Qp[setup.NS+k] + a*n[0][k])*dlambda2[i] + (dQp_dQ[setup.NS+k][i] - da_dQ[i]*n[0][k])*lambda3 + (Qp[setup.NS+k] - a*n[0][k])*dlambda3[i]);
		}
	}


	/*
	 * C. Total Energy
	 */
#pragma omp parallel for
	for (int i = 0; i <= setup.VAR-1; i++)
	{
		Jacob_inviscid[ne][i] = (dQp_dQ[ne][i]*A + Qp[ne]*dA_dQ[i] - dp_dQ[i]) * lambda1;
		Jacob_inviscid[ne][i] += (Qp[ne]*A - W[ne]) * dlambda1[i];
		Jacob_inviscid[ne][i] += dB_dQ[i] * ((Qp[ne]+Up[0]*a)*lambda2 + (Qp[ne]-Up[0]*a)*lambda3);
		Jacob_inviscid[ne][i] +=	B * ((dQp_dQ[ne][i] + dUa_dQ[i])*lambda2	+ (Qp[ne] + Up[0]*a)*dlambda2[i]);
		Jacob_inviscid[ne][i] +=	B * ((dQp_dQ[ne][i] - dUa_dQ[i])*lambda3	+ (Qp[ne] - Up[0]*a)*dlambda3[i]);
	}


	/*
	 * D. Other energy modes
	 */
#pragma omp parallel for
	for (int m = 1; m <= setup.NE-1; m++)
	{
		for (int i = 0; i <= setup.VAR-1; i++)
		{
			Jacob_inviscid[ne+m][i] =	(A*dQp_dQ[ne+m][i] + dA_dQ[i]*Qp[ne+m]) * lambda1;
			Jacob_inviscid[ne+m][i] +=	A*Qp[ne+m] * dlambda1[i];
			Jacob_inviscid[ne+m][i] += (B*dQp_dQ[ne+m][i] + dB_dQ[i]*Qp[ne+m]) * (lambda2 + lambda3);
			Jacob_inviscid[ne+m][i] += B*Qp[ne+m] * (dlambda2[i] + dlambda3[i]);
		}
	}
}
