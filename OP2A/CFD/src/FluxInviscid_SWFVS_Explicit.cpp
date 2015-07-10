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
#include <time.h>

#include "CFD/include/FluxInviscid.hpp"
#include "Math/include/OP2A_Matrix.hpp"
#include "Math/include/MathMisc.hpp"
#include "Common/include/Time_Info.hpp"
#include "Common/include/MultiDimension.hpp"

namespace OP2A{
namespace CFD{


void FluxInviscid::SWFVS_Explicit(Data::DataStorageVector<Data::DataStorage>& data1D_L, Data::DataStorageVector<Data::DataStorage>& data1D_R, CHEM::SpeciesSet& species_set, int ND,
									unsigned int type, unsigned int indexQ, unsigned int indexV, unsigned int indexW,
									vector< vector<double> >& normal_vector, int faceID,
									double dp, double dist_wall, double n_dot_wall, double alpha, double x0, double eps0,
									Data::DataStorage& Fn_inv)
{
	int VAR	= data1D_L(indexQ).numData;
	int CFD_NT	= omp_get_num_threads();
	int index_u	= species_set.NS;
	int index_E = species_set.NS + ND;
	int NE	= VAR - index_E;

	// Initialize Jacobians (Plus/minus)
	vector <Math::MATRIX> A_pm(2, Math::MATRIX(VAR, VAR, false));

	// 1. Caculate Omega (Pressure Switch factor)
	double omega	= 0.5 / (pow(alpha*dp, 2.0) + 1.0);
	double o_m_omega = 1.0 - omega;


	//  2. Calculate Flux at face (MAIN)
	for (int l = 0; l <= 1; l++)	// [0]: Plus, [1]:Minus
	{
		double signal	= pow(-1.0, l);

		Data::DataStorageVector<Data::DataStorage>	data1D_j(6);
		data1D_j(0).resize(VAR);	// Q
		data1D_j(1).resize(VAR);	// V
		data1D_j(2).resize(VAR);	// W
		data1D_j(3).resize(6);	// MIX
		data1D_j(4).resize(species_set.NS);	// Xs
		data1D_j(5).resize(species_set.NS);	// Ys

		// 2.1 Get Q at face Qj(+/-)
		switch(l)
		{
		case 0:
#pragma ivdep
			for(int i = 0; i <= VAR-1; i++)	data1D_j.data[0].data[i]	= o_m_omega*data1D_L(indexQ)(i)	+ omega*data1D_R(indexQ)(i);
			break;

		case 1:
#pragma ivdep
			for(int i = 0; i <= VAR-1; i++)	data1D_j.data[0].data[i]	= o_m_omega*data1D_R(indexQ)(i)	+ omega*data1D_L(indexQ)(i);
			break;
		}

		VariableChange::From_Q(type, data1D_j.data[0],  data1D_j.data[1],  data1D_j.data[2],  data1D_j.data[3], data1D_j.data[4],  data1D_j.data[5], species_set, ND, CFD_NT);


		// 2.2 Calculate additional variables
		double rho = data1D_j.data[3].data[0];

		double U_square = 0.0;
		for (int k = index_u; k <= index_u+ND-1; k++)	U_square += pow(data1D_j.data[1].data[k], 2.0);

		double H;	// Enthalpy
		H	= (data1D_j.data[0].data[index_E] + data1D_j.data[2].data[index_E]) / rho;

		vector<double> Up(ND, 0.0);

		for (int k = 0; k <= ND-1; k++)
		{
			for (int k1 = 0; k1 <= ND-1; k1++)
			{
				Up[k]	+= data1D_j.data[1].data[index_u+k1] * normal_vector[k][k1];
			}
		}


		// 2.3 Calculate dT/dQ, dp/dQ, a2
		Data::DataStorage	dp("dp_dQ", VAR);
		Data::DataStorageVector<Data::DataStorage>	dT(NE);
		double a2;
		double two_a2;
		double a;
		double two_a;

		for (int i = 0; i <= NE-1; i++)	dT(i).resize(VAR);
		Derivatives::dTdQ(data1D_j, species_set, ND, type, 0, 1, 2,	3, dT);
		Derivatives::dpdQ(data1D_j, dT, species_set, ND, type, 0, 1, 2, 3, dp);
		a2		= Derivatives::a2(data1D_j, dp, species_set, ND, type, 0, 1, 2, 3);

		two_a2	= 2.0 *a2;
		a		= sqrt(a2);
		two_a	= 2.0*a;


		// 2.4. Eigenvalues
		double eps;

		if (dist_wall < x0)	eps	= eps0 * n_dot_wall * (sqrt(U_square) + a);
		else				eps	= eps0 * (sqrt(U_square) + a);


		double  eps2	= eps * eps;
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

#pragma omp parallel for private(temp1, temp2)
		for (int s1 = 0; s1 <= species_set.NS-1; s1++)
		{
			// Species
			temp1	= -data1D_j.data[0].data[s1] / rho_a2 * 0.5;
			temp2	= data1D_j.data[5].data[s1] * Up[0] / two_a;

			for (int s2 = 0; s2 <= species_set.NS-1; s2++)
			{
				A_pm[l](s1, s2)	= temp1*dp.data[s2]*aux1 - temp2*aux2;
			}
			A_pm[l](s1, s1) += lambda1;


			// Momentum
			temp2	= data1D_j.data[5].data[s1] / two_a;

			for (int k2 = 0; k2 <= ND-1; k2++)
			{
				A_pm[l](s1, index_u + k2) = temp1*dp.data[index_u+k2]*aux1 + temp2*normal_vector[0][k2]*aux2;
			}


			// Energy modes
			for (int e2 = 0; e2 <= NE-1; e2++)
			{
				A_pm[l](s1, index_E + e2) = temp1*dp.data[index_E + e2]*aux1;
			}
		}


		// 2.5.2. Momentum
		for (int k1 = 0; k1 <= ND-1; k1++)
		{
			// Species
			temp1	= data1D_j.data[1].data[index_u + k1] / a2;
			temp2	= data1D_j.data[1].data[index_u + k1] * Up[0];

#pragma ivdep
			for (int s2 = 0; s2 <= species_set.NS-1; s2++)
			{
				A_pm[l](index_u+k1, s2)	= -(dp.data[s2]*temp1 - Up[0]*normal_vector[0][k1])*0.5*aux1
											+ (dp.data[s2]*normal_vector[0][k1] - temp2)/two_a*aux2;
			}


			// Momentum
			for (int k2 = 0; k2 <= ND-1; k2++)
			{
				A_pm[l](index_u+k1, index_u+k2)	= -(dp(index_u+k2)*temp1 + normal_vector[0][k1]*normal_vector[0][k2])*0.5*aux1
														+ (dp(index_u+k2)*normal_vector[0][k1] + data1D_j(1)(index_u+k1)*normal_vector[0][k2])/two_a*aux2;
			}
			A_pm[l](index_u+k1, index_u+k1)	+= lambda1;

			// Energy
			for (int e2 = 0; e2 <= NE-1; e2++)
			{
				A_pm[l](index_u+k1, index_E+e2) = -0.5*temp1*dp(index_E+e2)*aux1 + dp(index_E+e2)*normal_vector[0][k1]/two_a*aux2;
			}
		}


		// 2.5.3. Total Energy
		// Species
		temp1	= H/a2;
#pragma ivdep
		for (int s2 = 0; s2 <= species_set.NS-1; s2++)
		{
			A_pm[l](index_E, s2) = (-dp(s2)*temp1 + Up[0]*Up[0])*0.5*aux1	+ Up[0]/two_a*(dp(s2) - H)*aux2;
		}

		// Momentum
		for (int k2 = 0; k2 <= ND-1; k2++)
		{
			A_pm[l](index_E, index_u+k2) = -(dp(index_u+k2)*temp1 + Up[0]*normal_vector[0][k2])*0.5*aux1
												+ (dp(index_u+k2)*Up[0] + H*normal_vector[0][k2])/two_a*aux2;
		}

		// Energy
		for (int e2 = 0; e2 <= NE-1; e2++)
		{
			A_pm[l](index_E, index_E+e2) = -dp(index_E+e2)*temp1*0.5*aux1	+ dp(index_E+e2)*Up[0]/two_a*aux2;
		}
		A_pm[l](index_E, index_E) += lambda1;



		// 2.5.4. Other energy modes
		for (int e1 = 1; e1 <= NE-1; e1++)
		{
			// Species
			temp1	= data1D_j(0)(index_E+e1) / rho_a2;
			temp2	= data1D_j(0)(index_E+e1) / rho;

#pragma ivdep
			for (int s2 = 0; s2 <= species_set.NS-1; s2++)
			{
				A_pm[l](index_E+e1, s2)	= -dp(s2)*temp1*0.5*aux1 - Up[0]/two_a*temp2*aux2;
			}

			// Momentum
#pragma ivdep
			for (int k2 = 0; k2 <= ND-1; k2++)
			{
				A_pm[l](index_E+e1, index_u+k2) = -dp(index_u+k2)*temp1*0.5*aux1 + temp2*normal_vector[0][k2]/(2.0*a)*aux2;
			}

			// Energy
#pragma ivdep
			for (int e2 = 0; e2 <= NE-1; e2++)
			{
				A_pm[l](index_E+e1, index_E+e2) = -dp(index_E+e2)*temp1*0.5*aux1;
			}
			A_pm[l](index_E+e1, index_E+e1)	+= lambda1;
		}
	}

	 // 3. Calculate Flux at Face
#pragma omp parallel for
	for (int i = 0; i <= VAR-1; i++)
	{
		Fn_inv(i) = 0.0;
		for (int j = 0; j <= VAR-1; j++)
		{
			//Fn_inv(i)	+= A_pm[0](i,j)*data1D_L(indexQ)(j) + A_pm[1](i,j)*data1D_R(indexQ)(j);
			Fn_inv(i)	+= A_pm[0](i,j) * data1D_L.data[indexQ].data[j] + A_pm[1](i,j)*data1D_R.data[indexQ].data[j];
		}
	}
}

}
}