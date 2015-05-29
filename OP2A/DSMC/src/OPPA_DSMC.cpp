/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: Dec 15, 2014
 *      			Author: Minkwan Kim
 *
 * OPPA_DSMC.cpp
 * 			-  
 *  
 */


/*
 * ===============================
 * 		Include Header files
 * ===============================
 */
#include <omp.h>
#include <signal.h>
#include <ctime>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "../../Grid/include/OPPA_grid.hpp"
#include "../../Data_storage/include/OPPA_data.hpp"


#include "../include/DSMC_fns.hpp"

#include "../../utilities/include/general_fns.hpp"
#include "../../utilities/include/version.hpp"




using namespace std;




/*
 * =================================
 * 	Set Global variables
 * ================================
 */
int		*WHEREIS_CELL;					// Location of cells (cellptr/taskid)
int		NP					= 1;		// Total number of processor
int		P					= 0;		// Current processor number
int		OLD_OBJ				= 0;
int		NEW_OBJ				= 1;
int		FLAG_DIM			= 2;		// Flag for dimension (1 = Axisymmetric / 2 = 2D / 3 = 3D)
bool	FLAG_SAMPLING		= false;    // Flag for sampling
bool	FLAG_TERMINATION	= false;	// Flag for termination signal

double 	T_EPSIOLON;						// Relocation time tolerance

particle_object	*unused_object	= NULL;	// Unused object
#ifdef PARALLEL
	comm_packet	*unused_packet	= NULL; // Unused packages
#endif

Cell_type		**cell_ptr;				// Cell pointers array

int		n_sample			= 0; 		// Current number of sample step
int		n_object			= 0;		// Total number of object in simulation






int main( int argc, char *argv[] )
{
	/*
	 * ==============================================
	 * STEP 1: Initialize Parallel Communication:
	 * 		- Development Status : Done
	 *		- Last modified on: July 23, 2014
	 *						by: Minkwan Kim
	 * =============================================
	 * */
	time_t	t = time(0);  					// get time now
	struct	tm * now = localtime( & t );

	int		N_thread	= 24;
	double 	t0			= 0.0;
	double 	tstart;								// Simulation time
	double 	tstop;
	double 	t_simulation;

	tstart	= dtime();

#ifdef PARALLEL

#ifdef SP2
  mpc_environ(&numtask,&taskid);   /* Get number of parallel tasks */
#endif

#ifdef MPI
	MPI_Init(&argc,&argv);					// INITIALIZE MPI
	MPI_Comm_size(MPI_COMM_WORLD, &NP); 	// FINDOUT HOW MANY PROCESSORS IN THERE
	MPI_Comm_rank(MPI_COMM_WORLD, &P);  	// FINDOUT WHICH PROCESSOR I AM
	t0 = MPI_Wtime();
#endif

	if (NP > MAX_N_TASKS)	program_error_type1("Number of processors exceeds MAX_PROCESSOR. Need to adjust value of MAX_N_TASKS");
#endif



	/*
	 * =================================================
	 * STEP 2: Show Version information:
	 *		- Development Status : Done
	 *		- Last modified on: July 23, 2014
	 *		-				by: Minkwan Kim
	 * =================================================
	 */
	Ver_type ver;							// Version information class
	ver.primary		= 0;
	ver.secondary	= 0;
	ver.type		= "DSMC";
	ver.year		= now->tm_year + 1900;
	ver.month		= now->tm_mon + 1;
	ver.date		= now->tm_mday;

	if (P == 0)	ver.out_info();				// Show the version information





	/*
	 * =======================================================
	 * STEP 3: Read Problem information
	 * 		- Development Status: Version 1.0
	 * 		- Last modified on: July 23, 2014
	 * 						by: Minkwan Kim
	 * =======================================================
	 */
	OPPA_DSMC_setup problem;
	string mesh_file_name	= "RamC_2D_v1.cas";
	string output_file_name	= "test.plt";
	string output_title		= "Test_CASE";






	/*
	 * ======================================================================
	 * STEP 5: GRID GENERATION and/or READ (Unstructured Catersian grid)
	 * 		- Development Status: Version 1.0a
	 * 		- Last modified on: July 23, 2014
	 * 						by: Minkwan Kim
	 * ======================================================================
	 */
	GRID_CLASS grid;
	preprocessing_grid__read_process_ver1(mesh_file_name, FLUENT, false, 1000.0, grid);


	SOL_CLASS_BASIC	Cell_solution_Q;
	Cell_solution_Q.allocate_size(grid.NCM, 5);













	/*
	 * ===================================================================
	 * STEP 8: Preparing Output data
	 * 		- Development Status: Version 1.0
	 *====================================================================
	 */
		#ifdef MPI
			if (P == 0)
				cout << "  Write initial Output file  [t = " << MPI_Wtime()-t0 << endl;
		#else
			cout << "  Write initial Output file" << endl;
		#endif

		// Step 8.2 Write Initial data
		int	time_strand_ID = 0;
		data_print(P, grid, Cell_solution_Q, output_title, output_file_name, 0, time_strand_ID, TECPLOT, PRINT_NORMAL);

		#ifdef MPI
			if (P == 0)
				cout << "  --->Done  [t = " << MPI_Wtime()-t0 << endl;
		#else
			cout << "  --->Done" << endl;
		#endif








	/*
	 * =======================================================
	 * Physical Solving Loop (DSMC)
	 * 		- Development Status: Version 1.0
	 * 		- Last modified on: Dec/15/2014
	 * 						by: Minkwan Kim
	 * =======================================================
	 */
	double perf, perf2;

	if (P	== 0)	cout << "Start OPPA-DSMC...." << endl;
	while (problem.n_current < problem.n_total && FLAG_TERMINATION != true)
	{

		// 1. Reset Event counters
		problem.reset_counters();
		problem.n_current++;

		FLAG_SAMPLING	= false;
		if (problem.n_current > problem.n_startup	&&	problem.n_current % problem.intvl_sample == 0)
		{
			n_sample++;
			FLAG_SAMPLING	= true;
		}


		// 2. Calculate Cell
		calculate_cell(problem, cell_ptr, problem.grid.ND, perf, perf2);


	}






}



