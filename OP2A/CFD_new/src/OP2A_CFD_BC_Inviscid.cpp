/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 15, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_BC_Inviscid.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_CFD_BC.hpp"


int OP2A_BC_type_adjust(int BC_index)
{
	int bc_update;
	switch(BC_index)
	{
	case BC_WALL:
		bc_update	= 1;
		break;

	case BC_SYMMETRY:
		bc_update	= 1;
		break;

	case BC_AXIS:
		bc_update	= 1;
		break;

	case BC_ANODE:
		bc_update	= 1;
		break;

	case BC_CATHODE:
		bc_update	= 1;
		break;

	case BC_DIELECTRIC_WALL:
		bc_update	= 1;
		break;

	case BC_INLET:
		bc_update	= 2;
		break;

	case BC_FREESTREAM:
		bc_update	= 4;
		break;

	case BC_OUTLET:
		bc_update	= 3;
		break;
	}

	return (bc_update);
}



void CFD_assign_BC_inviscid(CFD_variable_setup_ver2 setup, vector < vector <double> > &Q_inlet, GRID_CLASS &grid,
							vector < vector<double> > &Qc, vector < vector<double> > &Qgc, int nt)
{
	int f, cr, cl;
	int f_ptr, cr_ptr, cl_ptr;

	kmp_set_defaults("KMP_AFFINITY = scatter");


#pragma omp parallel for private(f, cl, cr, f_ptr, cr_ptr, cl_ptr) num_threads(nt)
	for (int i = 0; i <= grid.NGM-1; i++)
	{
		f		= grid.cells_ghost.data_ptr[i]->face[0];
		f_ptr	= grid.faces.whereis[f];

		cl	= grid.faces.data_ptr[f_ptr]->cl[0];
		cr	= grid.faces.data_ptr[f_ptr]->cr[0];

		cl_ptr	= grid.cells.whereis[cl];
		cr_ptr	= grid.cells_ghost.whereis[-cr];


		int BC_index	= OP2A_BC_type_adjust(grid.cells_ghost.data_ptr[i]->BC);

		switch(BC_index)
		{
		case 1:
			CFD_BC_Inviscid_wall(setup, *(grid.cells.data_ptr[cl_ptr]), *(grid.faces.data_ptr[f_ptr]), *(grid.cells_ghost.data_ptr[cr_ptr]), Qc[cl_ptr], Qgc[cr_ptr]);
			break;

		case 2:
			CFD_BC_Inviscid_inlet(setup, *(grid.cells.data_ptr[cl_ptr]), *(grid.cells_ghost.data_ptr[cr_ptr]), Qc[cl_ptr], Qgc[cr_ptr], Q_inlet[1]);
			break;

		case 3:
			CFD_BC_Inviscid_exit(setup, *(grid.cells.data_ptr[cl_ptr]), *(grid.cells_ghost.data_ptr[cr_ptr]), Qc[cl_ptr], Qgc[cr_ptr]);
			break;

		case 4:
			CFD_BC_Inviscid_inlet(setup, *(grid.cells.data_ptr[cl_ptr]), *(grid.cells_ghost.data_ptr[cr_ptr]), Qc[cl_ptr], Qgc[cr_ptr], Q_inlet[0]);
			break;

		}
	}
}



void CFD_assign_BC_inviscid_complete(SOL_CFD &Solution_data, GRID_CLASS &grid, vector < vector <double> > &Q_IC, vector < vector <double> > &V_IC, SPECIES_DATA_BASIC &species, int nt)
{
	CFD_assign_BC_inviscid(Solution_data.setup, Q_IC, grid, Solution_data.Qc, Solution_data.Qgc,nt);

	// 2. Update V and Additional values of ghost cells
#pragma omp parallel for num_threads(nt)
	for (int i = 0; i <= grid.NGM-1; i++)
	{
		Solution_data.mixture_data_gc[i].calculate_data(Solution_data.setup.NS, Solution_data.setup.ID_T[ROT], Solution_data.Qgc[i], species.Cv_t, species.Cv_r, species.R, species.M);
		CFD_Q_to_V(Solution_data.Qgc[i], Solution_data.Vgc[i], Solution_data.setup, Solution_data.mixture_data_gc[i], species);
		CFD_V_to_W(Solution_data.Vgc[i], Solution_data.Wgc[i], Solution_data.setup, Solution_data.mixture_data_gc[i], species);
	}
}



void CFD_assign_BC_inviscid_implicit(CFD_variable_setup_ver2 setup, GRID_CLASS &grid, vector< vector< vector<double> > > &Jacobian_plus, vector< vector< vector<double> > > &Jacobian_minus, int nt)
{
	int f;
	int f_ptr;

//#pragma omp parallel for private(f, f_ptr) num_threads(nt)
	for (int i = 0; i <= grid.NGM-1; i++)
	{
		f				= grid.cells_ghost.data_ptr[i]->face[0];
		f_ptr			= grid.faces.whereis[f];
		int BC_index	= OP2A_BC_type_adjust(grid.cells_ghost.data_ptr[i]->BC);

		vector < vector<double>	> E			= vector_2D(setup.VAR, setup.VAR, 0.0);
		vector < vector<double>	> R			= vector_2D(setup.ND, setup.ND, 0.0);
		vector < vector<double>	> R_inv		= vector_2D(setup.ND, setup.ND, 0.0);
		for (int j = 0; j <= setup.VAR-1; j++)	E[j][j]	= 1.0;

		switch(BC_index)
		{
		case 1:
			for (int j = 0; j <= setup.ND-1; j++)
			{
				for (int k = 0; k <= setup.ND-1; k++)
				{
					R[j][k]		= grid.faces.data_ptr[f_ptr]->n[j][k];
					R_inv[j][k]	= grid.faces.data_ptr[f_ptr]->n[k][j];
				}
			}
			for (int k = 0; k <= setup.ND-1; k++)	R[0][k]	= -R[0][k];

			for (int k = 0; k <= setup.ND-1; k++)
			{
				for (int l = 0; l <= setup.ND-1; l++)
				{
					double aux = 0.0;
					for (int m = 0; m <= setup.ND-1; m++)
					{
						aux	+= R_inv[k][m] * R[m][l];
					}

					E[setup.NS + k][setup.NS + l]	= aux;
				}
			}

			for (int k = 0; k <= setup.VAR-1; k++)
			{
				for (int l = 0; l <= setup.VAR-1; l++)
				{
					double aux = 0.0;
					for (int m = 0; m <= setup.VAR-1; m++)
					{
						aux	+= Jacobian_minus[f_ptr][k][m] * E[m][l];
					}

					Jacobian_plus[f_ptr][k][l]	+= aux;
				}

				for (int l = 0; l <= setup.VAR-1; l++)		Jacobian_minus[f_ptr][k][l]	= 0.0;
			}
			break;

		case 2:
			for (int k = 0; k <= setup.VAR-1; k++)
			{
				for (int l = 0; l <= setup.VAR-1; l++)	Jacobian_minus[f_ptr][k][l]	= 0.0;
			}
			break;

		case 3:
			for (int k = 0; k <= setup.VAR-1; k++)
			{
				for (int l = 0; l <= setup.VAR-1; l++)	Jacobian_plus[f_ptr][k][l]	+= Jacobian_minus[f_ptr][k][l];
				for (int l = 0; l <= setup.VAR-1; l++)	Jacobian_minus[f_ptr][k][l]	= 0.0;
			}
			break;

		case 4:
			for (int k = 0; k <= setup.VAR-1; k++)
			{
				for (int l = 0; l <= setup.VAR-1; l++)	Jacobian_minus[f_ptr][k][l]	= 0.0;
			}
			break;
		}
	}
}


