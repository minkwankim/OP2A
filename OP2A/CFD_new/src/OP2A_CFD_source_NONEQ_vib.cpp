/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 28, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_source_NONEQ_vib.cpp
 * 			-  
 *  
 */




#include "../include/OP2A_CFD_Variable_change.hpp"
#include "../include/OP2A_CFD_reconstruct.hpp"
#include "../include/OP2A_CFD_Flux.hpp"


#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../../SETUP/include/OP2A_setup_problem.hpp"
#include "../../MATRIX/include/OPPA_Matrix.hpp"
#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../GRID/include/OP2A_grid.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"







void CFD_source_NONEQ_vib(GRID_CLASS &grid, SOL_CFD &Solution, SPECIES_DATA_BASIC &species, int nt)
{

	int ne	= Solution.setup.NS + Solution.setup.ND;

	for (int c = 0; c <= grid.NCM-1; c++)
	{
		double Q_tv	= 0.0;
		double Q_rv	= 0.0;
		double Q_ev = 0.0;
		double Q_v	= 0.0;

		double p 	= Solution.Wc[c][ne];
		double T	= Solution.Vc[c][ne];
		double Tv	= Solution.Vc[c][ne + Solution.setup.ID_T[VIB]];
		double Te	= Solution.Vc[c][ne + Solution.setup.ID_T[ELE]];

		double T_m_1_3	= pow(T, -1.0/3.0);

		for (int s = 0; s <= Solution.setup.NS-1; s++)
		{
			double tau_s = 0.0;

			int s_ID	= species.whereis[s];
			if((*species.data_entire)[s_ID].basic_data.type == MOLECULE)
			{
				double num	= 0.0;
				double den 	= 0.0;

				for (int r = 0; r <= Solution.setup.NS-1; r++)
				{
					int r_ID	= species.whereis[r];
					if((*species.data_entire)[r_ID].basic_data.type != ELECTRON)
					{
						double mu_sr	= species.M[s]*species.M[r] / (species.M[s] + species.M[r]);
						double theta_s	= min_vec((*species.data_entire)[s_ID].thermodynamic_properties.theta_vib, (*species.data_entire)[s_ID].thermodynamic_properties.n_vib_lvl);
						double A_sr		= 1.16e-3	* sqrt(mu_sr)	* pow(theta_s, 4.0/3.0);
						double B_sr		= 0.015 * pow(mu_sr, 0.25);
						double tau_sr	= 101325.0/p	* exp(A_sr*(T_m_1_3 - B_sr) - 18.42);

						num	+= Solution.mixture_data_c[c].Xs[r];
						den	+= Solution.mixture_data_c[c].Xs[r] / tau_sr;
					}
					else
					{
						double ne		= Solution.Qc[c][s]	/ (*species.data_entire)[r_ID].basic_data.m;
						if ((*species.data_entire)[s_ID].thermodynamic_properties.n_vib_ele > 0)
						{
							double tau_es 	= (*species.data_entire)[s_ID].tau_e_vib(Te, ne);
							Q_ev			= ((*species.data_entire)[r_ID].thermodynamic_properties.Cv[0]*Te - (*species.data_entire)[r_ID].thermodynamic_properties.Cv[0]*T) / tau_es;
						}
					}
				}

				tau_s = num/den;

				double sigma_s = 5.81E-21;
				double c_s = sqrt(8.*C_BOLTZMANN_SI*T*N_AVOGADRO_SI/PI/Solution.mixture_data_c[c].M_mix);
				double tau_park	= 1.0 / (sigma_s * c_s * (Solution.mixture_data_c[c].rho/Solution.mixture_data_c[c].M_mix*N_AVOGADRO_SI));

				tau_s	+= tau_park;

				double evs_star	= (*species.data_entire)[s_ID].e_vib(T);
				double evs		= (*species.data_entire)[s_ID].e_vib(Tv);

				Q_tv	+= Solution.Qc[c][s] * (evs_star - evs)/tau_s;
				Q_v		+= Solution.S_source[c][s] *(*species.data_entire)[s_ID].e_VEE(Tv);
			}
		}




		Solution.S_source[c][ne + Solution.setup.ID_T[VIB]]	+= Q_tv + Q_rv + Q_ev + Q_v;
	}
}
