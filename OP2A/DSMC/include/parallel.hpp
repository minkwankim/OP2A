/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: Dec 15, 2014
 *      			Author: Minkwan Kim
 *
 * parallel.hpp
 * 			-  
 *  
 */

#ifndef PARALLEL_HPP_
#define PARALLEL_HPP_

#ifdef PARALLEL

#include "particle.hpp"
#include "global_variables.hpp"

/*
 * Vendor Specific defines and includes
 */

#ifdef SP2
#include <mpproto.h>
#endif

#ifdef MPI
#include <mpi.h>
#endif

/*
 * Set parameters
 */
#define MAX_N_TASKS		1024		// Maximum number of task
#define MAX_N_MESG		100000		// Maximum number of object for MPI communication
#define MAX_N_PROC		24			// Maximum number of process for vector computing



/*
 * Message types
 */
#define OBJ_DATA     20
#define CELL_DATA   100
#define FIEL_DDATA  200
#define WALL_DATA   300
#define MESG_DATA   500




/*
 * Data type for message communication packet
 */
class comm_packet
{
public:
	  particle_type	prop;          			// Particle data

	  int    		old_cell_id;			// Old cell location of object
	  int			new_cell_id;  			// New cell location of object

	  double  		t_move;                	// Time left for movement
	  double		t_done;                	// Time done for movement

	  comm_packet 	*next;    				// Next object
};


extern comm_packet *unused_packet;  		// Holds unused packages
#endif

#endif /* PARALLEL_HPP_ */
