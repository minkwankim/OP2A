/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 15, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Flux_inviscid.cpp
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


void CFD_Flux_inviscid(GRID_CLASS &grid, SOL_CFD &Solution, SPECIES_DATA_BASIC &species, int order, int limiter, bool is_axis, double alpha, double x0, double eps0, int implicit, int nt)
{

	if (implicit != 0)
	{
#pragma omp parallel for num_threads(nt)
		for (int f = 0; f <= grid.NFM-1; f++)
		{
			int cl, cll, cr, crr;
			vector <double> x_f(Solution.setup.ND, 0.0);
			vector < vector <double> >	n = vector_2D(Solution.setup.ND, Solution.setup.ND, 0.0);

			for (int k = 0; k <= Solution.setup.ND-1; k++)	x_f[k]	= grid.faces.data_ptr[f]->x[k];
			for (int k1 = 0; k1 <= Solution.setup.ND-1; k1++)
				for (int k2 = 0; k2 <= Solution.setup.ND-1; k2++)
					n[k1][k2]	= grid.faces.data_ptr[f]->n[k1][k2];


			// Reconstruct Face left and right value
			vector<double>	W_cl(Solution.setup.VAR, 0.0);
			vector<double>	W_cll(Solution.setup.VAR, 0.0);
			vector<double>	W_cr(Solution.setup.VAR, 0.0);
			vector<double>	W_crr(Solution.setup.VAR, 0.0);

			vector<double>	x_cl(Solution.setup.ND, 0.0);
			vector<double>	x_cll(Solution.setup.ND, 0.0);
			vector<double>	x_cr(Solution.setup.ND, 0.0);
			vector<double>	x_crr(Solution.setup.ND, 0.0);


			cl	= 	grid.faces.data_ptr[f]->cl[0];
			cll	= 	grid.faces.data_ptr[f]->cl[1];

			cr	=	grid.faces.data_ptr[f]->cr[0];
			crr	=	grid.faces.data_ptr[f]->cr[1];


			// Right cell
			if (cr > 0)
			{
				int cr_ptr	= grid.cells.whereis[cr];

				W_cr	= Solution.Wc[cr_ptr];
				for (int k = 0; k <= Solution.setup.ND-1; k++)	x_cr[k]	= grid.cells.data_ptr[cr_ptr]->x[k];
			}
			else
			{
				int cr_ptr	= grid.cells_ghost.whereis[-cr];

				W_cr = Solution.Wgc[cr_ptr];
				for (int k = 0; k <= Solution.setup.ND-1; k++)	x_cr[k]	= grid.cells_ghost.data_ptr[cr_ptr]->x[k];
			}


			// Right-Right cell
			if (crr > 0)
			{
				int crr_ptr	= grid.cells.whereis[crr];

				W_crr	= Solution.Wc[crr_ptr];
				for (int k = 0; k <= Solution.setup.ND-1; k++)	x_crr[k]	= grid.cells.data_ptr[crr_ptr]->x[k];
			}
			else
			{
				int crr_ptr	= grid.cells_ghost.whereis[-crr];

				W_crr	= Solution.Wgc[crr_ptr];
				for (int k = 0; k <= Solution.setup.ND-1; k++)	x_crr[k]	= grid.cells_ghost.data_ptr[crr_ptr]->x[k];
			}


			// Left cell
			if (cl > 0)
			{
				int cl_ptr	= grid.cells.whereis[cl];

				W_cl	= Solution.Wc[cl_ptr];
				for (int k = 0; k <= Solution.setup.ND-1; k++)	x_cl[k]	= grid.cells.data_ptr[cl_ptr]->x[k];
			}
			else
			{
				int cl_ptr	= grid.cells_ghost.whereis[-cl];

				W_cl	= Solution.Wgc[cl_ptr];
				for (int k = 0; k <= Solution.setup.ND-1; k++)	x_cl[k]	= grid.cells_ghost.data_ptr[cl_ptr]->x[k];
			}

			// Left-left
			if (cll > 0)
			{
				int cll_ptr	= grid.cells.whereis[cll];

				W_cll	= Solution.Wc[cll_ptr];
				for (int k = 0; k <= Solution.setup.ND-1; k++)	x_cll[k] = grid.cells.data_ptr[cll_ptr]->x[k];
			}
			else
			{
				int cll_ptr	= grid.cells_ghost.whereis[-cll];

				W_cll	= Solution.Wgc[cll_ptr];
				for (int k = 0; k <= Solution.setup.ND-1; k++)	x_cll[k] = grid.cells_ghost.data_ptr[cll_ptr]->x[k];
			}

			// Reconstruct
			vector<double>	W_L(Solution.setup.VAR, 0.0);
			vector<double>	W_R(Solution.setup.VAR, 0.0);
			reconstruct_MUSCL_ver3(W_cll, W_cl, W_cr, W_crr, x_cll,	x_cl, x_cr, x_crr, x_f,	3.0/8.0, order, limiter, Solution.setup.ND, Solution.setup.VAR, W_L, W_R);

			vector<double>	V_L(Solution.setup.VAR, 0.0);
			vector<double>	V_R(Solution.setup.VAR, 0.0);

			vector<double>	Q_L(Solution.setup.VAR, 0.0);
			vector<double>	Q_R(Solution.setup.VAR, 0.0);

			CFD_mixture_data mixture_data_L;
			CFD_mixture_data mixture_data_R;

			mixture_data_L.calculate_data(Solution.setup.NS, Solution.setup.ID_T[ROT], W_L, species.Cv_t, species.Cv_r, species.R, species.M);
			mixture_data_R.calculate_data(Solution.setup.NS, Solution.setup.ID_T[ROT], W_R, species.Cv_t, species.Cv_r, species.R, species.M);

			CFD_W_to_V(W_L, V_L, Solution.setup, mixture_data_L, species);
			CFD_W_to_V(W_R, V_R, Solution.setup, mixture_data_R, species);

			CFD_V_to_Q(V_L, Q_L, Solution.setup, mixture_data_L, species);
			CFD_V_to_Q(V_R, Q_R, Solution.setup, mixture_data_R, species);

			double p_cl = W_cl[Solution.setup.NS + Solution.setup.ND];
			double p_cr = W_cr[Solution.setup.NS + Solution.setup.ND];
			double dp	= fabs(p_cl - p_cr)  /fmin(p_cl, p_cr);

			vector<double>	F_1_2(Solution.setup.VAR);

			double S = 0.0;
			if (is_axis	== true)	S	= grid.faces.data_ptr[f]->S * fabs(grid.faces.data_ptr[f]->x[1]);
			else					S	= grid.faces.data_ptr[f]->S;

			CFD_flux_1_2_FVM_single_ver3(Solution.setup, Q_L, V_L, mixture_data_L, Q_R, V_R, mixture_data_R, Solution.Flux_f_inviscid[f], n, S, species, Solution.Jacobian_inviscid_plus[f], Solution.Jacobian_inviscid_minus[f], dp, grid.faces.data_ptr[f]->dist_wall, grid.faces.data_ptr[f]->n_dot_wall, alpha, x0, eps0);
		}
	}
	else
	{
#pragma omp parallel for num_threads(nt)
		for (int f = 0; f <= grid.NFM-1; f++)
		{
			int cl, cll, cr, crr;
			vector <double> x_f(Solution.setup.ND, 0.0);
			vector < vector <double> >	n = vector_2D(Solution.setup.ND, Solution.setup.ND, 0.0);

			for (int k = 0; k <= Solution.setup.ND-1; k++)	x_f[k]	= grid.faces.data_ptr[f]->x[k];
			for (int k1 = 0; k1 <= Solution.setup.ND-1; k1++)
				for (int k2 = 0; k2 <= Solution.setup.ND-1; k2++)
					n[k1][k2]	= grid.faces.data_ptr[f]->n[k1][k2];


			// Reconstruct Face left and right value
			vector<double>	W_cl(Solution.setup.VAR, 0.0);
			vector<double>	W_cll(Solution.setup.VAR, 0.0);
			vector<double>	W_cr(Solution.setup.VAR, 0.0);
			vector<double>	W_crr(Solution.setup.VAR, 0.0);

			vector<double>	x_cl(Solution.setup.ND, 0.0);
			vector<double>	x_cll(Solution.setup.ND, 0.0);
			vector<double>	x_cr(Solution.setup.ND, 0.0);
			vector<double>	x_crr(Solution.setup.ND, 0.0);


			cl	= 	grid.faces.data_ptr[f]->cl[0];
			cll	= 	grid.faces.data_ptr[f]->cl[1];

			cr	=	grid.faces.data_ptr[f]->cr[0];
			crr	=	grid.faces.data_ptr[f]->cr[1];


			// Right cell
			if (cr > 0)
			{
				int cr_ptr	= grid.cells.whereis[cr];

				W_cr	= Solution.Wc[cr_ptr];
				for (int k = 0; k <= Solution.setup.ND-1; k++)	x_cr[k]	= grid.cells.data_ptr[cr_ptr]->x[k];
			}
			else
			{
				int cr_ptr	= grid.cells_ghost.whereis[-cr];

				W_cr = Solution.Wgc[cr_ptr];
				for (int k = 0; k <= Solution.setup.ND-1; k++)	x_cr[k]	= grid.cells_ghost.data_ptr[cr_ptr]->x[k];
			}


			// Right-Right cell
			if (crr > 0)
			{
				int crr_ptr	= grid.cells.whereis[crr];

				W_crr	= Solution.Wc[crr_ptr];
				for (int k = 0; k <= Solution.setup.ND-1; k++)	x_crr[k]	= grid.cells.data_ptr[crr_ptr]->x[k];
			}
			else
			{
				int crr_ptr	= grid.cells_ghost.whereis[-crr];

				W_crr	= Solution.Wgc[crr_ptr];
				for (int k = 0; k <= Solution.setup.ND-1; k++)	x_crr[k]	= grid.cells_ghost.data_ptr[crr_ptr]->x[k];
			}


			// Left cell
			if (cl > 0)
			{
				int cl_ptr	= grid.cells.whereis[cl];

				W_cl	= Solution.Wc[cl_ptr];
				for (int k = 0; k <= Solution.setup.ND-1; k++)	x_cl[k]	= grid.cells.data_ptr[cl_ptr]->x[k];
			}
			else
			{
				int cl_ptr	= grid.cells_ghost.whereis[-cl];

				W_cl	= Solution.Wgc[cl_ptr];
				for (int k = 0; k <= Solution.setup.ND-1; k++)	x_cl[k]	= grid.cells_ghost.data_ptr[cl_ptr]->x[k];
			}

			// Left-left
			if (cll > 0)
			{
				int cll_ptr	= grid.cells.whereis[cll];

				W_cll	= Solution.Wc[cll_ptr];
				for (int k = 0; k <= Solution.setup.ND-1; k++)	x_cll[k] = grid.cells.data_ptr[cll_ptr]->x[k];
			}
			else
			{
				int cll_ptr	= grid.cells_ghost.whereis[-cll];

				W_cll	= Solution.Wgc[cll_ptr];
				for (int k = 0; k <= Solution.setup.ND-1; k++)	x_cll[k] = grid.cells_ghost.data_ptr[cll_ptr]->x[k];
			}

			// Reconstruct
			vector<double>	W_L(Solution.setup.VAR, 0.0);
			vector<double>	W_R(Solution.setup.VAR, 0.0);
			reconstruct_MUSCL_ver3(W_cll, W_cl, W_cr, W_crr, x_cll,	x_cl, x_cr, x_crr, x_f,	3.0/8.0, order, limiter, Solution.setup.ND, Solution.setup.VAR, W_L, W_R);

			vector<double>	V_L(Solution.setup.VAR, 0.0);
			vector<double>	V_R(Solution.setup.VAR, 0.0);

			vector<double>	Q_L(Solution.setup.VAR, 0.0);
			vector<double>	Q_R(Solution.setup.VAR, 0.0);

			CFD_mixture_data mixture_data_L;
			CFD_mixture_data mixture_data_R;

			mixture_data_L.calculate_data(Solution.setup.NS, Solution.setup.ID_T[ROT], W_L, species.Cv_t, species.Cv_r, species.R, species.M);
			mixture_data_R.calculate_data(Solution.setup.NS, Solution.setup.ID_T[ROT], W_R, species.Cv_t, species.Cv_r, species.R, species.M);

			CFD_W_to_V(W_L, V_L, Solution.setup, mixture_data_L, species);
			CFD_W_to_V(W_R, V_R, Solution.setup, mixture_data_R, species);

			CFD_V_to_Q(V_L, Q_L, Solution.setup, mixture_data_L, species);
			CFD_V_to_Q(V_R, Q_R, Solution.setup, mixture_data_R, species);

			double p_cl = W_cl[Solution.setup.NS + Solution.setup.ND];
			double p_cr = W_cr[Solution.setup.NS + Solution.setup.ND];
			double dp	= fabs(p_cl - p_cr)  /fmin(p_cl, p_cr);

			vector<double>	F_1_2(Solution.setup.VAR);

			double S = 0.0;
			if (is_axis	== true)	S	= grid.faces.data_ptr[f]->S * fabs(grid.faces.data_ptr[f]->x[1]);
			else					S	= grid.faces.data_ptr[f]->S;

			CFD_flux_1_2_FVM_single_ver2(Solution.setup, Q_L, V_L, mixture_data_L, Q_R, V_R, mixture_data_R, Solution.Flux_f_inviscid[f], n, S, species, dp, grid.faces.data_ptr[f]->dist_wall, grid.faces.data_ptr[f]->n_dot_wall, alpha, x0, eps0);
			//CFD_flux_1_2_AUSM_single(Solution.setup, Q_L, V_L, W_L, mixture_data_L, Q_R, V_R, W_R, mixture_data_R, Solution.Flux_f_inviscid[f], n, S, species);

		}
	}
}
