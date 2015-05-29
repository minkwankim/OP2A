/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 15, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_FLux_Inviscid_FVM.cpp
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


void CFD_flux_1_2_FVM_single_ver3(CFD_variable_setup_ver2 setup,
								vector<double>	&Q_L, vector<double> &U_L, CFD_mixture_data &variable_data_L,
								vector<double>	&Q_R, vector<double> &U_R, CFD_mixture_data &variable_data_R,
								vector<double>	&F_1_2, vector< vector<double> > &n, double S, SPECIES_DATA_BASIC	&species,
								vector< vector<double> >	&Jacobian_plus, vector< vector<double> >	&Jacobian_minus,
								double dp, double dist_wall, double n_dot_wall,
								double alpha, double x0, double eps0)
{

	// Initialize variables: Jacobians
	vector <vector <vector <double> > > A_pm(2);
	A_pm[0]	= vector_2D(setup.VAR, setup.VAR, 0.0);
	A_pm[1]	= vector_2D(setup.VAR, setup.VAR, 0.0);


	for (int i = 0; i <= setup.VAR-1; i++)
	{
		for (int j = 0; j <= setup.VAR-1; j++)
		{
			Jacobian_plus[i][j]		= 0.0;
			Jacobian_minus[i][j]	= 0.0;
		}
	}


	// 1. Calculate omega(pressure switch)
	double omega		= 0.5 / ((alpha*dp)*(alpha*dp) + 1.0);
	double one_m_omega	= 1.0 - omega;


	/*
	 * 2. Calculate Flux Part (MAIN)
	 */
	for (int l	= 0; l <= 1; l++)	// [0]:plus	[1]: minux
	{
		//	+: POSITIVE JACOBIAN   -: NEGATIVE JACOBIAN
		double signal = pow(-1.0, l);


		// 1. Calculate Qj and Uj from Q_L and Q_R
		vector<double>		Qj(setup.VAR, 0.0);
		vector<double>		Vj(setup.VAR, 0.0);
		vector<double>		Wj(setup.VAR, 0.0);
		CFD_mixture_data 	variable_data_j;

		// Calculate Qj and Uj
		switch(l)
		{
		case 0:
			for (int i = 0; i <= setup.VAR-1; i++)	Qj[i]	= one_m_omega*Q_L[i]	+ omega*Q_R[i];
			break;
		case 1:
			for (int i = 0; i <= setup.VAR-1; i++)	Qj[i]	= one_m_omega*Q_R[i]	+ omega*Q_L[i];
			break;
		}

		variable_data_j.calculate_data(setup.NS, setup.ID_T[ROT], Qj, species.Cv_t, species.Cv_r, species.R, species.M);
		CFD_Q_to_V(Qj, Vj, setup, variable_data_j, species);
		CFD_V_to_W(Vj, Wj, setup, variable_data_j, species);


		// Additional variables
		double U_square	= 0.0;
		for (int k = 0; k <= setup.ND-1; k++)	U_square += pow(Vj[setup.NS+k], 2.0);

		double	H	= (Qj[setup.NS+setup.ND] + Wj[setup.NS+setup.ND]) / variable_data_j.rho;


		// Calculate Normal/Tangential/Tangential2 velocities
		vector<double>	Up(setup.ND, 0.0);

		for (int k1 = 0; k1 <= setup.ND-1; k1++)
		{
			for (int k = 0; k <= setup.ND-1; k++)
			{
				Up[k1]	+= Vj[setup.NS+k] * n[k1][k];
			}
		}



		// 2. Calculate dT_dQ
		vector < vector <double> >	dT_dQ = vector_2D(setup.NE, setup.VAR, 0.0);
		CFD_Calculate_dT_dQ(setup, Qj, Vj, variable_data_j, dT_dQ, species);


		// 3. Calculate dp_dQ
		vector<double>	dp_dQ(setup.VAR, 0.0);
		CFD_Calculate_dp_dQ(setup, Vj, variable_data_j, dT_dQ, dp_dQ,	species);

		double a2		= CFD_Calculate_a2(setup, Qj, Wj, variable_data_j, dp_dQ, species);
		double two_a2	= 2.0 *a2;

		double a;
		if (a2 >= 0.0)	a	= sqrt(a2);
		else			program_error_type1("Error negative a2: Please check dp and dT calculation!!");



		// 4. Eigenvalues
		double eps;
		if (dist_wall < x0)	eps	= eps0 * n_dot_wall * (fabs(Up[0]) + a);
		else				eps	= eps0 * (fabs(Up[0]) + a);

		double 	lambda1	= 0.5 * (Up[0] 		+ signal*sqrt(Up[0]*Up[0] 		+ eps*eps));
		double 	lambda2	= 0.5 * ((Up[0]+a)	+ signal*sqrt(pow(Up[0]+a, 2.0) + eps*eps));
		double	lambda3	= 0.5 * ((Up[0]-a)	+ signal*sqrt(pow(Up[0]-a, 2.0) + eps*eps));


		// 5. Calculate A+-
		double aux1	= 2.0*lambda1 - lambda2 - lambda3;
		double aux2	= lambda2 - lambda3;
		double temp1;
		double temp2;

		// 1. Species
		for (int s = 0; s <= setup.NS-1; s++)
		{
			// Species
			temp1	= -Qj[s] / (variable_data_j.rho * a2) * 0.5;
			temp2	= (Qj[s] / variable_data_j.rho) * Up[0] / (2.0*a);
			for (int i = 0; i <= setup.NS-1; i++)	A_pm[l][s][i]	= temp1*dp_dQ[i]*aux1 - temp2*aux2;
			A_pm[l][s][s]	+= lambda1;

			// Momentum
			temp2	= (Qj[s] / variable_data_j.rho) / (2.0*a);
			for (int i = 0; i <= setup.ND-1; i++)	A_pm[l][s][setup.NS+i]	= temp1*dp_dQ[setup.NS+i]*aux1 + temp2*n[0][i]*aux2;

			// Energy modes
			for (int i = 0; i <= setup.NE-1; i++)	A_pm[l][s][setup.NS+setup.ND+i]	= temp1*dp_dQ[setup.NS+setup.ND+i]*aux1;
		}


		// 2. Momentum
		for (int k = 0; k <= setup.ND-1; k++)
		{
			// Species
			temp1	= Vj[setup.NS+k]/a2;
			temp2	= Vj[setup.NS+k] * Up[0];
			for (int i = 0; i <= setup.NS-1; i++)	A_pm[l][setup.NS+k][i]	= -(dp_dQ[i]*temp1 - Up[0]*n[0][k])*0.5*aux1 + (dp_dQ[i]*n[0][k] - temp2)/(2.0*a)*aux2;

			// Momentum
			for (int i = 0; i <= setup.ND-1; i++)	A_pm[l][setup.NS+k][setup.NS+i]	= -(dp_dQ[setup.NS+i]*temp1 + n[0][k]*n[0][i])*0.5*aux1 + (dp_dQ[setup.NS+i]*n[0][k] + Vj[setup.NS+k]*n[0][i])/(2.0*a)*aux2;
			A_pm[l][setup.NS+k][setup.NS+k]	+= lambda1;

			// Energy
			for (int i = 0; i <= setup.NE-1; i++)	A_pm[l][setup.NS+k][setup.NS+setup.ND+i] = -0.5*temp1*dp_dQ[setup.NS+setup.ND+i]*aux1 + dp_dQ[setup.NS+setup.ND+i]*n[0][k]/(2.0*a)*aux2;
		}


		// 3. Total Energy
		// Species
		temp1	= H/a2;
		for (int i = 0; i <= setup.NS-1; i++)	A_pm[l][setup.NS+setup.ND][i]	= (-dp_dQ[i]*temp1 + Up[0]*Up[0])*0.5*aux1	+ Up[0]/(2.0*a)*(dp_dQ[i] - H)*aux2;

		// Momentum
		for (int i = 0; i <= setup.ND-1; i++)	A_pm[l][setup.NS+setup.ND][setup.NS+i]	= -(dp_dQ[setup.NS+i]*temp1 + Up[0]*n[0][i])*0.5*aux1	+ (dp_dQ[setup.NS+i]*Up[0] + H*n[0][i])/(2.0*a)*aux2;

		// Energy
		for (int i = 0; i <= setup.NE-1; i++)	A_pm[l][setup.NS+setup.ND][setup.NS+setup.ND+i]	= -dp_dQ[setup.NS+setup.ND+i]*temp1*0.5*aux1	+ dp_dQ[setup.NS+setup.ND+i]*Up[0]/(2.0*a)*aux2;
		A_pm[l][setup.NS+setup.ND][setup.NS+setup.ND]	+= lambda1;

		// 4. Other energy modes
		for (int k = 1; k <= setup.NE-1; k++)
		{
			// Species
			temp1	= Qj[setup.NS+setup.ND+k]/variable_data_j.rho/a2;
			temp2	= Qj[setup.NS+setup.ND+k]/variable_data_j.rho;
			for (int i = 0; i <= setup.NS-1; i++)	A_pm[l][setup.NS+setup.ND+k][i]	= -dp_dQ[i]*temp1*0.5*aux1	- Up[0]/(2.0*a)*temp2*aux2;

			// Momentum
			for (int i = 0; i <= setup.ND-1; i++)	A_pm[l][setup.NS+setup.ND+k][setup.NS+i] = -dp_dQ[setup.NS+i]*temp1*0.5*aux1	+ temp2*n[0][i]/(2.0*a)*aux2;

			// Energy
			for (int i = 0; i <= setup.NE-1; i++)	A_pm[l][setup.NS+setup.ND+k][setup.NS+setup.ND+i] = -dp_dQ[setup.NS+setup.ND+i]*temp1*0.5*aux1;
			A_pm[l][setup.NS+setup.ND+k][setup.NS+setup.ND+k]	+= lambda1;
		}



		switch(l)
		{
		case 0:
			CFD_Jacobian_inviscid_single(setup,	Qj, Vj, Wj, variable_data_j, dT_dQ, dp_dQ, Jacobian_plus, signal, n, species, dp, dist_wall, n_dot_wall, x0, eps0);
			for (int i = 0; i <= setup.VAR-1; i++)
				for (int j = 0; j <= setup.VAR-1; j++)
					Jacobian_plus[i][j]	= Jacobian_plus[i][j]*S;
			break;
		case 1:
			CFD_Jacobian_inviscid_single(setup,	Qj, Vj, Wj, variable_data_j, dT_dQ, dp_dQ, Jacobian_minus, signal, n, species, dp, dist_wall, n_dot_wall, x0, eps0);
			for (int i = 0; i <= setup.VAR-1; i++)
				for (int j = 0; j <= setup.VAR-1; j++)
					Jacobian_minus[i][j]	= Jacobian_minus[i][j]*S;
			break;
		}
	}



	/*
	 * ========================================================
	 * [Final] Calculate Flux at face
	 * ========================================================
	 */
	for (int i = 0; i <= setup.VAR-1; i++)	F_1_2[i]	= 0.0;		// Initialize it

	for (int i = 0; i <= setup.VAR-1; i++)
	{
		for (int j = 0; j <= setup.VAR-1; j++)
		{
			F_1_2[i]	+= A_pm[0][i][j]*Q_L[j] + A_pm[1][i][j]*Q_R[j];
		}
	}
}





void CFD_flux_1_2_FVM_single_ver2(CFD_variable_setup_ver2 setup,
								vector<double>	&Q_L, vector<double> &U_L, CFD_mixture_data &variable_data_L,
								vector<double>	&Q_R, vector<double> &U_R, CFD_mixture_data &variable_data_R,
								vector<double>	&F_1_2, vector< vector<double> > &n, double S, SPECIES_DATA_BASIC	&species,
								double dp, double dist_wall, double n_dot_wall,
								double alpha, double x0, double eps0)
{

	// Initialize variables: Jacobians
	vector <vector <vector <double> > > A_pm(2);
	A_pm[0]	= vector_2D(setup.VAR, setup.VAR, 0.0);
	A_pm[1]	= vector_2D(setup.VAR, setup.VAR, 0.0);


	// 1. Calculate omega(pressure switch)
	double omega		= 0.5 / ((alpha*dp)*(alpha*dp) + 1.0);
	double one_m_omega	= 1.0 - omega;


	/*
	 * 2. Calculate Flux Part (MAIN)
	 */
	for (int l	= 0; l <= 1; l++)	// [0]:plus	[1]: minux
	{
		//	+: POSITIVE JACOBIAN   -: NEGATIVE JACOBIAN
		double signal = pow(-1.0, l);


		// 1. Calculate Qj and Uj from Q_L and Q_R
		vector<double>		Qj(setup.VAR, 0.0);
		vector<double>		Vj(setup.VAR, 0.0);
		vector<double>		Wj(setup.VAR, 0.0);
		CFD_mixture_data 	variable_data_j;

		// Calculate Qj and Uj
		switch(l)
		{
		case 0:
			for (int i = 0; i <= setup.VAR-1; i++)	Qj[i]	= one_m_omega*Q_L[i]	+ omega*Q_R[i];
			break;
		case 1:
			for (int i = 0; i <= setup.VAR-1; i++)	Qj[i]	= one_m_omega*Q_R[i]	+ omega*Q_L[i];
			break;
		}

		variable_data_j.calculate_data(setup.NS, setup.ID_T[ROT], Qj, species.Cv_t, species.Cv_r, species.R, species.M);
		CFD_Q_to_V(Qj, Vj, setup, variable_data_j, species);
		CFD_V_to_W(Vj, Wj, setup, variable_data_j, species);


		// Additional variables
		double U_square	= 0.0;
		for (int k = 0; k <= setup.ND-1; k++)	U_square += pow(Vj[setup.NS+k], 2.0);

		double	H	= (Qj[setup.NS+setup.ND] + Wj[setup.NS+setup.ND]) / variable_data_j.rho;


		// Calculate Normal/Tangential/Tangential2 velocities
		vector<double>	Up(setup.ND, 0.0);

		for (int k1 = 0; k1 <= setup.ND-1; k1++)
		{
			for (int k = 0; k <= setup.ND-1; k++)
			{
				Up[k1]	+= Vj[setup.NS+k] * n[k1][k];
			}
		}



		// 2. Calculate dT_dQ
		vector < vector <double> >	dT_dQ = vector_2D(setup.NE, setup.VAR, 0.0);
		CFD_Calculate_dT_dQ(setup, Qj, Vj, variable_data_j, dT_dQ, species);


		// 3. Calculate dp_dQ
		vector<double>	dp_dQ(setup.VAR, 0.0);
		CFD_Calculate_dp_dQ(setup, Vj, variable_data_j, dT_dQ, dp_dQ,	species);

		double a2		= CFD_Calculate_a2(setup, Qj, Wj, variable_data_j, dp_dQ, species);
		double two_a2	= 2.0 *a2;

		double a;
		if (a2 >= 0.0)	a	= sqrt(a2);
		else			program_error_type1("Error negative a2: Please check dp and dT calculation!!");


		// 4. Eigenvalues
		double eps;
		if (dist_wall < x0)	eps	= eps0 * n_dot_wall * (fabs(Up[0]) + a);
		else				eps	= eps0 * (fabs(Up[0]) + a);

		double 	lambda1	= 0.5 * (Up[0] 		+ signal*sqrt(Up[0]*Up[0] 		+ eps*eps));
		double 	lambda2	= 0.5 * ((Up[0]+a)	+ signal*sqrt(pow(Up[0]+a, 2.0) + eps*eps));
		double	lambda3	= 0.5 * ((Up[0]-a)	+ signal*sqrt(pow(Up[0]-a, 2.0) + eps*eps));




		// 5. Calculate A+-
		double aux1	= 2.0*lambda1 - lambda2 - lambda3;
		double aux2	= lambda2 - lambda3;
		double temp1;
		double temp2;

		// 1. Species
		for (int s = 0; s <= setup.NS-1; s++)
		{
			// Species
			temp1	= -Qj[s] / (variable_data_j.rho * a2) * 0.5;
			temp2	= (Qj[s] / variable_data_j.rho) * Up[0] / (2.0*a);
			for (int i = 0; i <= setup.NS-1; i++)	A_pm[l][s][i]	= temp1*dp_dQ[i]*aux1 - temp2*aux2;
			A_pm[l][s][s]	+= lambda1;

			// Momentum
			temp2	= (Qj[s] / variable_data_j.rho) / (2.0*a);
			for (int i = 0; i <= setup.ND-1; i++)	A_pm[l][s][setup.NS+i]	= temp1*dp_dQ[setup.NS+i]*aux1 + temp2*n[0][i]*aux2;

			// Energy modes
			for (int i = 0; i <= setup.NE-1; i++)	A_pm[l][s][setup.NS+setup.ND+i]	= temp1*dp_dQ[setup.NS+setup.ND+i]*aux1;
		}


		// 2. Momentum
		for (int k = 0; k <= setup.ND-1; k++)
		{
			// Species
			temp1	= Vj[setup.NS+k]/a2;
			temp2	= Vj[setup.NS+k] * Up[0];
			for (int i = 0; i <= setup.NS-1; i++)	A_pm[l][setup.NS+k][i]	= -(dp_dQ[i]*temp1 - Up[0]*n[0][k])*0.5*aux1 + (dp_dQ[i]*n[0][k] - temp2)/(2.0*a)*aux2;

			// Momentum
			for (int i = 0; i <= setup.ND-1; i++)	A_pm[l][setup.NS+k][setup.NS+i]	= -(dp_dQ[setup.NS+i]*temp1 + n[0][k]*n[0][i])*0.5*aux1 + (dp_dQ[setup.NS+i]*n[0][k] + Vj[setup.NS+k]*n[0][i])/(2.0*a)*aux2;
			A_pm[l][setup.NS+k][setup.NS+k]	+= lambda1;

			// Energy
			for (int i = 0; i <= setup.NE-1; i++)	A_pm[l][setup.NS+k][setup.NS+setup.ND+i] = -0.5*temp1*dp_dQ[setup.NS+setup.ND+i]*aux1 + dp_dQ[setup.NS+setup.ND+i]*n[0][k]/(2.0*a)*aux2;
		}


		// 3. Total Energy
		// Species
		temp1	= H/a2;
		for (int i = 0; i <= setup.NS-1; i++)	A_pm[l][setup.NS+setup.ND][i]	= (-dp_dQ[i]*temp1 + Up[0]*Up[0])*0.5*aux1	+ Up[0]/(2.0*a)*(dp_dQ[i] - H)*aux2;

		// Momentum
		for (int i = 0; i <= setup.ND-1; i++)	A_pm[l][setup.NS+setup.ND][setup.NS+i]	= -(dp_dQ[setup.NS+i]*temp1 + Up[0]*n[0][i])*0.5*aux1	+ (dp_dQ[setup.NS+i]*Up[0] + H*n[0][i])/(2.0*a)*aux2;

		// Energy
		for (int i = 0; i <= setup.NE-1; i++)	A_pm[l][setup.NS+setup.ND][setup.NS+setup.ND+i]	= -dp_dQ[setup.NS+setup.ND+i]*temp1*0.5*aux1	+ dp_dQ[setup.NS+setup.ND+i]*Up[0]/(2.0*a)*aux2;
		A_pm[l][setup.NS+setup.ND][setup.NS+setup.ND]	+= lambda1;

		// 4. Other energy modes
		for (int k = 1; k <= setup.NE-1; k++)
		{
			// Species
			temp1	= Qj[setup.NS+setup.ND+k]/variable_data_j.rho/a2;
			temp2	= Qj[setup.NS+setup.ND+k]/variable_data_j.rho;
			for (int i = 0; i <= setup.NS-1; i++)	A_pm[l][setup.NS+setup.ND+k][i]	= -dp_dQ[i]*temp1*0.5*aux1	- Up[0]/(2.0*a)*temp2*aux2;

			// Momentum
			for (int i = 0; i <= setup.ND-1; i++)	A_pm[l][setup.NS+setup.ND+k][setup.NS+i] = -dp_dQ[setup.NS+i]*temp1*0.5*aux1	+ temp2*n[0][i]/(2.0*a)*aux2;

			// Energy
			for (int i = 0; i <= setup.NE-1; i++)	A_pm[l][setup.NS+setup.ND+k][setup.NS+setup.ND+i] = -dp_dQ[setup.NS+setup.ND+i]*temp1*0.5*aux1;
			A_pm[l][setup.NS+setup.ND+k][setup.NS+setup.ND+k]	+= lambda1;
		}
	}



	/*
	 * ========================================================
	 * [Final] Calculate Flux at face
	 * ========================================================
	 */
	for (int i = 0; i <= setup.VAR-1; i++)	F_1_2[i]	= 0.0;		// Initialize it

	for (int i = 0; i <= setup.VAR-1; i++)
	{
		for (int j = 0; j <= setup.VAR-1; j++)
		{
			F_1_2[i]	+= A_pm[0][i][j]*Q_L[j] + A_pm[1][i][j]*Q_R[j];
		}
	}
}


