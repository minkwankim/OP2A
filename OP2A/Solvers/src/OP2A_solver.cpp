/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 11, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_solver.cpp
 * 			-  
 *  
 */


#include <iostream>
#include <fstream>

#include "Common/include/Array2D.hpp"

#include "Common/include/Common.hpp"
#include "Common/include/Time_StopWatch.hpp"


#include "Configure/include/Config.hpp"
#include "Configure/include/ConfigArgs.hpp"






using namespace std;
using namespace OP2A;
using namespace OP2A::Common;
using namespace OP2A::Config;




int main(int argc, char** argv)
{
	CPUTime time_running;
	time_running.initStartTime();




	/*
	 * ==============================================
	 * STEP 1: Initialize Parallel Communication:
	 * 		- Development Status : Done
	 *		- Last modified on: July 23, 2014
	 *						by: Minkwan Kim
	 * =============================================
	 * */
	/*
	time_t	t = time(0);  					// get time now
	struct	tm * now = localtime( & t );

	double 	t0			= 0.0;
	double 	tstart;								// Simulation time
	double 	tstop;
	double 	t_simulation	= 0.0;
	double 	dt				= 0.0;

	tstart	= dtime();
	*/
#ifdef PARALLEL

#ifdef SP2
	mpc_environ(&numtask, &taskid);   /* Get number of parallel tasks */
#endif

#ifdef MPI
	MPI_Init(&argc,&argv);					// INITIALIZE MPI
	MPI_Comm_size(MPI_COMM_WORLD, &NP); 	// FINDOUT HOW MANY PROCESSORS IN THERE
	MPI_Comm_rank(MPI_COMM_WORLD, &P);  	// FINDOUT WHICH PROCESSOR I AM
	t0 = MPI_Wtime();
#endif

	if (NP > MAX_N_TASKS)	program_error_type1("Number of processors exceeds MAX_PROCESSOR. Need to adjust value of MAX_N_TASKS");
#endif
	//kmp_set_defaults("KMP_AFFINITY = scatter");



	//cout << temp2 <<endl;
	time_running.takeStopTime();
	double t0 = time_running.getDeltaT();



	cout << "HAHAHHA END" << endl;
	return (0);
}
