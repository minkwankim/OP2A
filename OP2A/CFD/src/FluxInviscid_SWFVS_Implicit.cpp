/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 26, 2015
 *      			Author: Minkwan Kim
 *
 * FluxInviscid_SWFVS_Implicit.cpp
 * 			-  
 *  
 */


#include <omp.h>

#include "CFD/include/FluxInviscid.hpp"
#include "Math/include/OP2A_Matrix.hpp"
#include "Math/include/MathMisc.hpp"


namespace OP2A{
namespace CFD{


void FluxInviscid::SWFVS_Implicit(Data::DataStorageVector<Data::DataStorage>& data1D_L, Data::DataStorageVector<Data::DataStorage>& data1D_R, CHEM::SpeciesSet& species_set, int ND,
									unsigned int type, unsigned int indexQ, unsigned int indexV, unsigned int indexW,
									vector< vector<double> >& normal_vector,
									double dp, double dist_wall, double n_dot_wall, double alpha, double x0, double eps0,
									Data::DataStorage& Fn_inv)
{
	int VAR	= data1D_L(indexQ).numData;
	int CFD_NT	= omp_get_num_threads();
	int index_u	= species_set.NS;
	int index_E = species_set.NS + ND;
	int NE	= VAR - index_E;

	// Initialize Jacobians (Plus/minus)
	vector<Math::MATRIX> A_pm(2, Math::MATRIX(VAR, VAR, false));


	// 1. Caculate Omega (Pressure Switch factor)
	double omega	= 0.5 / (pow(alpha*dp, 2.0) + 1.0);
	double o_m_omega = 1.0 - omega;


	//  2. Calculate Flux at face (MAIN)
	for (int l = 0; l <= 1; l++)	// [0]: Plus, [1]:Minus
	{
		double signal	= pow(-1.0, l);

		Data::DataStorageVector<Data::DataStorage>	data1D_j(3);
		data1D_j(0).resize(VAR);	// Q
		data1D_j(1).resize(VAR);	// V
		data1D_j(2).resize(VAR);	// W


		// 2.1 Get Q at face Qj(+/-)
		switch(l)
		{
		case 0:
#pragma omp parallel for
			for(int i = 0; i <= VAR-1; i++)	data1D_j(0)(i)	= o_m_omega*data1D_L(indexQ)(i)	+ omega*data1D_R(indexQ)(i);
			break;

		case 1:
#pragma omp parallel for
			for(int i = 0; i <= VAR-1; i++)	data1D_j(0)(i)	= o_m_omega*data1D_R(indexQ)(i)	+ omega*data1D_L(indexQ)(i);
			break;
		}

		VariableChange::Q_to_V(type, CFD_NT, data1D_j(0), species_set, ND, data1D_j(1));
		VariableChange::V_to_W(type, CFD_NT, data1D_j(1), species_set, ND, data1D_j(2));


		// 2.2 Calculate additional variables
		double rho = 0.0;
#pragma omp parallel for reduction(+: rho)
		for (int s = 0; s <= species_set.NS-1; s++)	rho += data1D_j(0)(s);

		double U_square = 0.0;
		for (int k = index_u; k <= index_u+ND-1; k++)	U_square += pow(data1D_j(1)(k), 2.0);

		double H;	// Enthalpy
		H	= (data1D_j(0)(index_E) + data1D_j(2)(index_E)) / rho;

		vector<double> Up(ND, 0.0);

		for (int k = 0; k <= ND-1; k++)
		{
			for (int k1 = 0; k1 <= ND-1; k1++)
			{
				Up[k]	+= data1D_j(1)(index_u+k1) * normal_vector[k][k1];
			}
		}

		// 2.3 Calculate dT/dQ, dp/dQ, a2
		Data::DataStorage	dp("dp_dQ", VAR);
		Data::DataStorageVector<Data::DataStorage>	dT(NE);
		for (int i = 0; i <= NE-1; i++)	dT(i).resize(VAR);

		Derivatives::dTdQ(data1D_j, species_set, ND, type, 0, 1, 2,	dT);
		Derivatives::dpdQ(data1D_j, dT, species_set, ND, type, 0, 1, 2, dp);

		double a2;
		double two_a2;
		double a;

		a2		= Derivatives::a2(data1D_j, dp, species_set, ND, type, 0, 1, 2);
		two_a2	= 2.0 *a2;
		a		= sqrt(a2);


		// 2.4. Eigenvalues
		double eps;

		if (dist_wall < x0)	eps	= eps0 * n_dot_wall * (Math::fabs<double>(Up[0]) + a);
		else				eps	= eps0 * (Math::fabs<double>(Up[0]) + a);

		double eps2	= eps * eps;
		double 	lambda1	= 0.5 * (Up[0] 		+ signal*sqrt(Up[0]*Up[0] 		+ eps2));
		double 	lambda2	= 0.5 * ((Up[0]+a)	+ signal*sqrt(pow(Up[0]+a, 2.0) + eps2));
		double	lambda3	= 0.5 * ((Up[0]-a)	+ signal*sqrt(pow(Up[0]-a, 2.0) + eps2));


		// 2.5. Calculate A (+/-)
		double aux1	= 2.0*lambda1 - lambda2 - lambda3;
		double aux2	= lambda2 - lambda3;
		double temp1;
		double temp2;

		// 2.5.1 Species
		double rho_a2	= rho * a2;

		for (int s = 0; s <= species_set.NS-1; s++)
		{
			// Species
			temp1	= -data1D_j(0)(s) / rho_a2 * 0.5;
			temp2	= data1D_j(0)(s)/rho * Up[0] / (2.0*a);

			for (int index_s = 0; index_s <= species_set.NS-1; index_s++)
			{
				A_pm[l](s, index_s)	= temp1*dp(index_s)*aux1 - temp2*aux2;
			}
			A_pm[l](s,s) += lambda1;


			// Momentum
			temp2	= data1D_j(0)(s)/rho / (2.0*a);
			for (int index_k = 0; index_k <= ND-1; index_k++)
			{
				A_pm[l](s, index_u + index_k) = temp1*dp(index_u+index_k)*aux1 + temp2*normal_vector[0][index_k]*aux2;
			}


			// Energy modes
			for (int index_e = 0; index_e <= NE-1; index_e++)
			{
				A_pm[l](s, index_E + index_e) = temp1*dp(index_E + index_e)*aux1;
			}
		}


		// 2.5.2. Momentum
		for (int k = 0; k <= ND-1; k++)
		{
			// Species
			temp1	= data1D_j(1)(index_u + k) / a2;
			temp2	= data1D_j(1)(index_u + k) * Up[0];

			for (int index_s = 0; index_s <= species_set.NS-1; index_s++)
			{
				A_pm[l](index_u+k, index_s)	= -(dp(index_s)*temp1 - Up[0]*normal_vector[0][k])*0.5*aux1 + (dp(index_s)*normal_vector[0][k] - temp2)/(2.0*a)*aux2;
			}


			// Momentum
			for (int index_k = 0; index_k <= ND-1; index_k++)
			{
				A_pm[l](index_u+k, index_u+index_k)	= -(dp(index_u+k)*temp1 + normal_vector[0][k]*normal_vector[0][index_k])*0.5*aux1
														+ (dp(index_u+index_k)*normal_vector[0][k] + data1D_j(1)(index_u+k)*normal_vector[0][index_k])/(2.0*a)*aux2;
			}
			A_pm[l](index_u+k, index_u+k)	+= lambda1;

			// Energy
			for (int index_e = 0; index_e <= NE-1; index_e++)
			{
				A_pm[l](index_u+k, index_E+index_e) = -0.5*temp1*dp(index_E+index_e)*aux1 + dp(index_E+index_e)*normal_vector[0][k]/(2.0*a)*aux2;
			}
		}


		// 2.5.3. Total Energy
		// Species
		temp1	= H/a2;
		for (int index_s = 0; index_s <= species_set.NS-1; index_s++)
		{
			A_pm[l](index_E, index_s) = (-dp(index_s)*temp1 + Up[0]*Up[0])*0.5*aux1	+ Up[0]/(2.0*a)*(dp(index_s) - H)*aux2;
		}

		// Momentum
		for (int index_k = 0; index_k <= ND-1; index_k++)
		{
			A_pm[l](index_E, index_u+index_k) = -(dp(index_u+index_k)*temp1 + Up[0]*normal_vector[0][index_k])*0.5*aux1
												+ (dp(index_u+index_k)*Up[0] + H*normal_vector[0][index_k])/(2.0*a)*aux2;
		}

		// Energy
		for (int index_e = 0; index_e <= NE-1; index_e++)
		{
			A_pm[l](index_E, index_E+index_e) = -dp(index_E+index_e)*temp1*0.5*aux1	+ dp(index_E+index_e)*Up[0]/(2.0*a)*aux2;
		}
		A_pm[l](index_E, index_E) += lambda1;



		// 2.5.4. Other energy modes
		for (int e = 1; e <= NE-1; e++)
		{
			// Species
			temp1	= data1D_j(0)(index_E+e) / rho_a2;
			temp2	= data1D_j(0)(index_E+e) / rho;

			for (int index_s = 0; index_s <= species_set.NS-1; index_s++)
			{
				A_pm[l](index_E+e, index_s)	= -dp(index_s)*temp1*0.5*aux1 - Up[0]/(2.0*a)*temp2*aux2;
			}

			// Momentum
			for (int index_k = 0; index_k <= ND-1; index_k++)
			{
				A_pm[l](index_E+e, index_u+index_k) = -dp(index_u+index_k)*temp1*0.5*aux1 + temp2*normal_vector[0][index_k]/(2.0*a)*aux2;
			}

			// Energy
			for (int index_e = 0; index_e <= NE-1; index_e++)
			{
				A_pm[l](index_E+e, index_E+index_e) = -dp(index_E+index_e)*temp1*0.5*aux1;
			}
			A_pm[l](index_E+e, index_E+e)	+= lambda1;
		}
	}



	 // 3. Calculate Flux at Face
	for (int i = 0; i <= VAR-1; i++)	Fn_inv(i) =	0.0;

	Fn_inv.data	= A_pm[0]*data1D_L(indexQ).data + A_pm[1]*data1D_R(indexQ).data;
}

}
}
