/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 7, 2015
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


void FluxInviscid::SWFVS_Implicit(Data::DataStorageVector<Data::DataStorage>& data1D_L, Data::DataStorageVector<Data::DataStorage>& data1D_R, CHEM::SpeciesSet& species_set, int ND,
									unsigned int type, unsigned int indexQ, unsigned int indexV, unsigned int indexW,
									vector< vector<double> >& normal_vector, int faceID,
									double dp, double dist_wall, double n_dot_wall, double alpha, double x0, double eps0,
									Data::DataStorage& Fn_inv,
									Data::DataStorage2D& dFdQ_plus,
									Data::DataStorage2D& dFdQ_minus)
{
	int VAR	= data1D_L(indexQ).numData;
	int CFD_NT	= omp_get_num_threads();
	int index_u	= species_set.NS;
	int index_E = species_set.NS + ND;
	int NE	= VAR - index_E;

	// Initialize Jacobians (Plus/minus)
	vector <Math::MATRIX> A_pm(2, Math::MATRIX(VAR, VAR, false));
	Data::DataStorage2D	dFdQpm("dFdQ", VAR, VAR);

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

//		if (dist_wall < x0)	eps	= eps0 * n_dot_wall * (sqrt(U_square) + a);
//		else				eps	= eps0 * (sqrt(U_square) + a);
		if (dist_wall < x0)	eps	= eps0 * n_dot_wall * (Math::fabs<double>(Up[0]) + a);
		else				eps	= eps0 * (Math::fabs<double>(Up[0]) + a);


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

#pragma ivdep
		for (int s1 = 0; s1 <= species_set.NS-1; s1++)
		{
			// Species
			double temp1A	= -data1D_j.data[0].data[s1] / rho_a2 * 0.5;
			double temp2A	= data1D_j.data[5].data[s1] * Up[0] / two_a;

			for (int s2 = 0; s2 <= species_set.NS-1; s2++)
			{
				A_pm[l](s1, s2)	= temp1A*dp.data[s2]*aux1 - temp2A*aux2;
			}
			A_pm[l](s1, s1) += lambda1;


			// Momentum
			temp2A	= data1D_j.data[5].data[s1] / two_a;

			for (int k2 = 0; k2 <= ND-1; k2++)
			{
				A_pm[l](s1, index_u + k2) = temp1A*dp.data[index_u+k2]*aux1 + temp2A*normal_vector[0][k2]*aux2;
			}


			// Energy modes
			for (int e2 = 0; e2 <= NE-1; e2++)
			{
				A_pm[l](s1, index_E + e2) = temp1A*dp.data[index_E + e2]*aux1;
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
#pragma ivdep
			for (int k2 = 0; k2 <= ND-1; k2++)
			{
				A_pm[l](index_u+k1, index_u+k2)	= -(dp(index_u+k2)*temp1 + normal_vector[0][k1]*normal_vector[0][k2])*0.5*aux1
														+ (dp(index_u+k2)*normal_vector[0][k1] + data1D_j(1)(index_u+k1)*normal_vector[0][k2])/two_a*aux2;
			}
			A_pm[l](index_u+k1, index_u+k1)	+= lambda1;

			// Energy
#pragma ivdep
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
#pragma ivdep
		for (int k2 = 0; k2 <= ND-1; k2++)
		{
			A_pm[l](index_E, index_u+k2) = -(dp(index_u+k2)*temp1 + Up[0]*normal_vector[0][k2])*0.5*aux1
												+ (dp(index_u+k2)*Up[0] + H*normal_vector[0][k2])/two_a*aux2;
		}

		// Energy
#pragma ivdep
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
				A_pm[l](index_E+e1, index_u+k2) = -dp(index_u+k2)*temp1*0.5*aux1 + temp2*normal_vector[0][k2]/two_a*aux2;
			}

			// Energy
#pragma ivdep
			for (int e2 = 0; e2 <= NE-1; e2++)
			{
				A_pm[l](index_E+e1, index_E+e2) = -dp(index_E+e2)*temp1*0.5*aux1;
			}
			A_pm[l](index_E+e1, index_E+e1)	+= lambda1;
		}



		/*
		 * Calculate
		 */
		// 1. Calculate d2p_dQ
		Data::DataStorage2D	d2p("dp_dQ", VAR, VAR);
		Derivatives::d2pdQ(data1D_j, dT, species_set, ND, type, 0, 1, 2, 3, d2p);

		// 2. Calculate da2_dQ and da_dQ
		Data::DataStorage	da2("da2", VAR);
		Data::DataStorage	da("da", VAR);
		Derivatives::da2dQ(data1D_j, dp, d2p, species_set, ND, 0, 1, 2, 3, da2);
		Derivatives::dadQ(da2, a, da);

		// 3. Calculate dU_dQ
		Data::DataStorage	dU_dQ("dU_dQ", VAR);
		double Up_rho = Up[0] / rho;

#pragma ivdep
		for (int s = 0; s <= species_set.NS-1; s++)
		{
			dU_dQ(s) = -Up_rho;
		}

#pragma ivdep
		for (int k = 0; k <= ND-1; k++)
		{
			dU_dQ(species_set.NS+k) = normal_vector[0][k]/rho;
		}

#pragma ivdep
		for (int e = species_set.NS+ND; e <= VAR-1; e++)
		{
			dU_dQ(e) = 0.0;
		}


		// 4. Calculate dlambda_dQ
		Data::DataStorage	dlambda1("dlambda1", VAR);
		Data::DataStorage	dlambda2("dlambda2", VAR);
		Data::DataStorage	dlambda3("dlambda3", VAR);

		double Up_Up2_p_eps2;
		double Uppa_Uppa2_p_eps2;
		double Upma_Upma2_p_eps2;

		if (Up[0] == 0 && eps2 == 0.0)
		{
			Up_Up2_p_eps2	= 1.0;
		}
		else
		{
			Up_Up2_p_eps2	= Up[0] / sqrt(Up[0]*Up[0] + eps2);
		}

		Uppa_Uppa2_p_eps2 	= (Up[0] + a) / sqrt((Up[0]+a)*(Up[0]+a) + eps2);
		Upma_Upma2_p_eps2	= (Up[0] - a) / sqrt((Up[0]-a)*(Up[0]-a) + eps2);

#pragma ivdep
		for (int i = 0; i <= VAR-1; i++)
		{
			dlambda1(i)	= 0.5 * (1.0 + signal*Up_Up2_p_eps2) * dU_dQ(i);
			dlambda2(i)	= 0.5 * (1.0 + signal*Uppa_Uppa2_p_eps2) * (dU_dQ(i) + da(i));
			dlambda3(i)	= 0.5 * (1.0 + signal*Upma_Upma2_p_eps2) * (dU_dQ(i) - da(i));
		}


		/*
		 * [PRE] Calculate predefined variables
		 */
		double p	= data1D_j(2)(index_E);
		double A	= p/a2 * dp(index_E);
		double B	= 0.5 * (rho - A);
		Data::DataStorage dA_dQ("dA_dQ", VAR);
		Data::DataStorage dB_dQ("dB_dQ", VAR);
		Data::DataStorage drho_dQ("drho_dQ", VAR);

		double dpdE_a4	= dp(index_E) /a2/a2;
		double p_a2		= p / a2;

#pragma ivdep
		for (int i = 0; i <= VAR-1;	i++)
		{
			dA_dQ(i)	= (dp(i)*a2 - da2(i)*p)*dpdE_a4	+ p_a2*d2p(index_E, i);
		}

#pragma ivdep
		for (int s = 0; s <= species_set.NS-1; 	s++)
		{
			drho_dQ(s)	= 1.0;
		}

#pragma ivdep
		for (int s = species_set.NS; s <= VAR-1; s++)
		{
			drho_dQ(s)	= 0.0;
		}

#pragma ivdep
		for (int i = 0; i <= VAR-1;	i++)
		{
			dB_dQ(i)	= 0.5 * (drho_dQ(i) - dA_dQ(i));
		}

		Data::DataStorage Qp("Qp", VAR);
#pragma ivdep
		for (int i = 0; i <= VAR-1;	i++)
		{
			Qp(i)	= data1D_j(0)(i) / rho;
		}
		Qp(index_E)	+= p/rho;


		Data::DataStorage2D dQp_dQ("dQp_dQ", VAR, VAR);
		double rho2	= rho * rho;
#pragma ivdep
		for (int r = 0; r <= VAR-1;	r++)
		{
#pragma ivdep
			for (int s = 0; s <= VAR-1;	s++)
			{
				dQp_dQ(r, s) = (Math::Delta_fn<double>(s,r)*rho	- drho_dQ(s)*data1D_j(0)(r)) / rho2;
			}
		}

#pragma ivdep
		for (int s = 0; s <= VAR-1;	s++)
		{
			dQp_dQ(index_E, s)	+= (dp(s)*rho	- drho_dQ(s)*p) / rho2;
		}

		Data::DataStorage dUa_dQ("dUa_dQ", VAR);
#pragma ivdep
		for (int i = 0; i <= VAR-1;	i++)
		{
			dUa_dQ(i)	= dU_dQ(i)*a	+ Up[0]*da(i);
		}


		/*
		 * A. Species
		 */
#pragma ivdep
		for (int s = 0; s <= species_set.NS-1; s++)
		{
			for (int i = 0; i <= VAR-1; i++)
			{
				dFdQpm(s,i) = (A*dQp_dQ(s,i) + dA_dQ(i)*Qp(s))*lambda1	+ A*Qp(s)*dlambda1(i)
							+ (B*dQp_dQ(s,i) + dB_dQ(i)*Qp(s))*(lambda2 + lambda3)	+ B*Qp(s)*(dlambda2(i) + dlambda3(i));
			}
		}

		/*
		 * B. Momentum
		 */

		for (int k = 0; k <= ND-1; k++)
		{
#pragma ivdep
			for (int i = 0; i <= VAR-1; i++)
			{
				dFdQpm(index_u+k,i)	 = (A*dQp_dQ(index_u+k,i) + dA_dQ(i)*Qp(index_u+k))*lambda1	+ A*Qp(index_u+k)*dlambda1(i);
				dFdQpm(index_u+k,i)	+= dB_dQ(i) * ((Qp(index_u+k) + a*normal_vector[0][k])*lambda2 + (Qp(index_u+k) - a*normal_vector[0][k])*lambda3);
				dFdQpm(index_u+k,i)	+= B * ((dQp_dQ(index_u+k,i) + da(i)*normal_vector[0][k])*lambda2 + (Qp(index_u+k) + a*normal_vector[0][k])*dlambda2(i) + (dQp_dQ(index_u+k,i) - da(i)*normal_vector[0][k])*lambda3 + (Qp(index_u+k) - a*normal_vector[0][k])*dlambda3(i));
			}
		}

		/*
		 * C. Total Energy
		 */
#pragma ivdep
		for (int i = 0; i <= VAR-1; i++)
		{
			dFdQpm(index_E,i)	 = (dQp_dQ(index_E, i)*A + Qp(index_E)*dA_dQ(i) - dp(i)) * lambda1;
			dFdQpm(index_E,i)	+= (Qp(index_E)*A - p) * dlambda1(i);
			dFdQpm(index_E,i)	+= dB_dQ(i) * ((Qp(index_E) + Up[0]*a)*lambda2 + (Qp(index_E) - Up[0]*a)*lambda3);
			dFdQpm(index_E,i)	+= B * ((dQp_dQ(index_E,i) + dUa_dQ(i))*lambda2	+ (Qp(index_E) + Up[0]*a)*dlambda2(i));
			dFdQpm(index_E,i)	+= B * ((dQp_dQ(index_E,i) - dUa_dQ(i))*lambda3	+ (Qp(index_E) - Up[0]*a)*dlambda3(i));
		}

		/*
		 * D. Other energy modes
		 */
		for (int e1 = index_E+1; e1 <= NE-1; e1++)
		{
#pragma ivdep
			for (int i = 0; i <= VAR-1; i++)
			{
				dFdQpm(e1,i)	 = (A*dQp_dQ(e1,i) + dA_dQ(i)*Qp(e1))*lambda1;
				dFdQpm(e1,i)	+= A*Qp(e1) * dlambda1(i);
				dFdQpm(e1,i)	+= (B*dQp_dQ(e1,i) + dB_dQ(i)*Qp(e1)) * (lambda2 + lambda3);
				dFdQpm(e1,i)	+= B*Qp(e1) * (dlambda2(i) + dlambda3(i));
			}
		}

		if (l == 0)
		{
			dFdQ_plus	= dFdQpm;
		}
		else
		{
			dFdQ_minus	= dFdQpm;
		}
	}

	 // 3. Calculate Flux at Face
#pragma ivdep
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
