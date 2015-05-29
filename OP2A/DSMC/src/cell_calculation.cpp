/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: Dec 15, 2014
 *      			Author: Minkwan Kim
 *
 * cell_calculation.cpp
 * 			-  
 *  
 */


#include "../include/DSMC_fns.hpp"
#include "../../utilities/include/general_fns.hpp"
#include "../../mesh/include/cell_DSMC.hpp"


void calculate_cell(OPPA_DSMC_setup &problem, Cell_type	**cell_ptr, int ND, double &perf, double &perf2)
{
	int				ncpu;
	Cell_type		*cell;
	Cell_type		*newcell;
	particle_object	*obj;
	particle_type	particles[MAX_NUM_OBJECT_PER_CELL];



	// 1. Initialize for Parallel(MPI)
#ifdef PARALLEL
	Cell_type	*cell_current;

	comm_packet	*packet;
	comm_packet	buffer[MAX_N_MESG]
	comm_packet	*object_buffer[MAX_N_TASKS];

	int	num_packet_send[MAX_N_TASKS];
	int	num_packet_recv[MAX_N_TASKS];

	for (int n = 0; n < problem.NP; n++)	num_packet_send[n]	= 0.0;

#ifdef MPI
	MPI_Request	request_send[MAX_N_TASKS];
	MPI_Request	request_recv[MAX_N_TASKS];
	MPI_Status	status;
#endif

#endif


	// 2. Reset variables
	clock_t t_calc	= 0.0;
	clock_t t_comm	= 0.0;
	clock_t	t_now	= clock();

	problem.reset_counters();



	// 3. Calculate Force Field (PIC)




	// For all cells
	int cell_ID, cell_ID_old, cell_ID_new;
	int num_object, num_object_old, num_object_new;

	for (int i = 0; i <= problem.NCM-1; i++)
	{
		cell		= cell_ptr[i];					// Current cell
		cell_ID		= cell->base.cell_ID;			// Cell ID
		num_object	= cell->base.num_object[OLD_OBJ];	// Number of objects(particles) in cell
		obj		= cell->base.object[OLD_OBJ];


		// 4.1 Copy data into Static array
		for (int j = 0; j <= num_object-1; j++)
		{
			particle_object *object_temp	= obj;

			particles[j]	= object_temp->prop;
			obj			= obj->next;
			ADD<particle_object>(object_temp, unused_object);
		}

		cell->base.num_object[OLD_OBJ]	= 0;
		cell->base.object[OLD_OBJ]		= NULL;


		// 4.2 Update Effective weight for cell and local time step
		double t_move;
		double t_done;
		double t_s	= problem.t_ref 		* problem.T_SCALE[cell_ID];
		double Wp	= problem.W_particle 	* problem.WEIGHT[cell_ID] 		* problem.T_SCALE[cell_ID];



		// 4.3A For Collision (Need to include)



		// 4.3B For Sampling (Need to include)


		// 4.3C For Chemistry (Need to include)




		// 4.4 Generate new particles in cells with inflow BC
		num_object_old	= num_object;
		num_object_new	= num_object;


		// 4.5 Move particles in cell
		cell_ID_old	= cell_ID;
		for (int j	= 0; j <= num_object_new-1; j++)
		{
			GET<particle_object>(obj, unused_object);
			obj->prop	= particles[j];

			if (j < num_object_old)	t_move	= t_s;
			else					t_move	= t_s *rand();

			t_done	= 0.0;

			cell_ID_new	= cell->base.cell_ID;
			move_particle(cell_ptr, obj, t_move, t_done, cell_ID_new, problem.M_move, problem.M_bad_move, ND);


			// Evaluate the movement
			if (cell_ID_new > 0)
			{
				newcell	= cell_ptr[WHEREIS_CELL[cell_ID_new]];

				// Clone or destroy particles


			}
			else
			{
				if (cell_ID_new == CELL_ID_REMOVE)
				{
					ADD<particle_object>(obj, unused_object);
				}
				else
				{
#ifdef	PARALLEL
					cell_ID_new	= -cell_ID_new;
					ncpu		= WHEREIS_CELL[cell_ID_new];
					GET<comm_packet>(packet, unused_packet);

					packet->particles	= obj->particles;
					packet->t_move		= t_move;
					packet->t_done		= t_done;
					packet->old_cell_id	= cell_ID_old;
					packet->new_cell_id	= cell_ID_new;

					ADD<comm_packet>(packet, unused_packet, object_buffer[ncpu]);
					ADD<particle_object>(obj, unused_object);
					num_packet_send[ncpu]++;
#endif
				}
			}
		}
	}

#ifdef	PARALLEL



#endif


	t_calc	+= clock() - t_now;
	perf	= 0.0;
	perf2	= 0.0;

	if (problem.n_particle_tot	> 0)
	{
		perf	= (t_calc + t_comm) / problem.n_particle_tot;

		if (perf != 0.0)
			perf2	= (t_calc/problem.n_particle_tot) / perf * 100;
	}

}



