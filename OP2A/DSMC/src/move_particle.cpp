/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: Dec 18, 2014
 *      			Author: Minkwan Kim
 *
 * move_particle.cpp
 * 			-  
 *  
 */


#include "../include/constants.hpp"
#include "../include/global_variables.hpp"
#include "../include/particle_object.hpp"

#include "../../utilities/include/general_fns.hpp"
#include "../../utilities/include/error_code.hpp"

#include "../../mesh/include/cell_DSMC.hpp"


#define ERROR_CODE_PRIMARY	ERROR_CODE_MOVE_PARTICLE



int find_face_to_heat(int ND, int num_face, double (*n)[3], double (*Xcf)[3], double *X, double *V, double &t_hit, int &f_hit)
{
	int find = 0;

	t_hit	= MATH_INFINITY;
	f_hit	= -1;

	for (int f = 0; f <= num_face-1; f++)
	{
		double n_dot_X	= 0.0;
		double n_dot_V	= 0.0;
		double t;

		for (int i = 0; i <= ND-1; i++)
		{
			n_dot_X	+= n[f][i] * (Xcf[f][i] - X[i]);
			n_dot_V	+= n[f][i] * V[i];
		}

		if (n_dot_V > 0.0)
		{
			t	= n_dot_X / n_dot_V;

			if (t < t_hit)
			{
				t_hit 	= t;
				f_hit	= f;
			}
		}
	}

	if (f_hit	== -1)	find = 1;

	return (find);
}






void move_particle(Cell_type **cell_p, particle_object *object, double &t_move, double &t_done, int &cell_ID_new, int &num_move, int &num_bad_move, int ND)
{
	int 				num_relocate	= 0;
	Cell_type			*cell;
	particle_type		*particle;
	Error_message_type	error_message;


	// 1. Assign start cell
	cell	= cell_p[WHEREIS_CELL[cell_ID_new]];


	// 2. Assign particle properties
	particle	= &(object->prop);

	double X[3];
	double V[3];
	for (int i = 0; i <= ND-1; i++)
	{
		X[i]	= particle->X[i];
		V[i]	= particle->V[i];
	}


	// 3. Movement
	int 	nface;
	int 	f_hit;
	int		is_find = 0;
	double 	t_hit;

	for (int l	= 0; l < MAX_NUM_EDGES; l++)
	{
		nface	= cell->geometry.num_face;


		while(is_find != 1 && num_relocate < 5)
		{
			is_find	= find_face_to_heat(ND, nface, cell->geometry.n, cell->geometry.Xcf, X, V, t_hit, f_hit);

			if (is_find != 1)
			{
				num_relocate++;
				for (int i = 0; i <= ND-1; i++)
				{
					X[i]	= 0.99*X[i]	+ 0.01*cell->geometry.xc[i];
				}
			}
		}


		// Error
		if (is_find != 1)
		{
			error_message.primary_code		= ERROR_CODE_MOVE_PARTICLE;;
			error_message.secondary_code	= ERROR_CODE_NO_FIND;

			error_message.location_primary		= cell_ID_new;
			error_message.location_primary_name	= "Cell ID";

			error_message.location_secondary		= 0;
			error_message.location_secondary_name	= "NONE";

			error_message.message	= "Check the position of particle";

			error_message.print_message();
		}



		// Case for movement within this cell
		if (t_hit > t_move)
		{
			for (int i = 0; i <= ND-1; i++)
			{
				particle->X[i]	= X[i] + V[i]*t_move;
			}

			cell_ID_new	= cell->base.cell_ID;
			return;
		}


		// Case: crossing cell boundary
		num_move++;
		if (cell->base.neighbor[f_hit] > 0)
		{
			t_move	= t_hit + (t_move - t_hit)*cell->base.t_ratio[f_hit];

			if (WHEREIS_CELL[cell->base.neighbor[f_hit]] >= 0)
			{
				cell	= cell_p[WHEREIS_CELL[cell->base.neighbor[f_hit]]];
			}
			else
			{
				cell_ID_new	= -cell->base.neighbor[f_hit] -1;
				return;
			}
		}
		else	// BOUNDARY CONDITION
		{
			int wall;
			int bc_flag = 0;

			switch(-cell->base.neighbor[f_hit])
			{
			case BC_INLET:
				bc_flag = 1;
				break;

			case BC_OUTLET:
				bc_flag = 1;
				break;

			case BC_FREESTREAM:
				bc_flag = 1;
				break;
			}


			if (bc_flag == 1)
			{
				cell_ID_new	= CELL_ID_REMOVE;
				return;
			}


			///
			// Need to include Wall BCs
			///
			for (int i = 0; i <= ND-1; i++)	particle->X[i]	= X[i] + V[i]*t_hit;


			for (int i = 0; i <= ND-1; i++)	V[i]	= particle->V[i];
			for (int i = 0; i <= ND-1; i++)	X[i]	= particle->X[i];
			t_move	= t_move - t_hit;
		}

	}



	// Bad movement
	num_bad_move++;
}









