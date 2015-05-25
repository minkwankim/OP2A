/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 25, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_Application.cpp
 * 			-  
 *  
 */


#include "Common/include/CodeLocation.hpp"
#include "Common/include/Exception_NPExceed.hpp"


#include "../include/OP2A_Application.hpp"


/*
 * FN-01: Prepararization OP2A
 * @author	Minkwan Kim
 * @version 1.0	25/5/2015
 */
void ApplicationOP2A::preparation(int argc, char *argv[], string modulename)
{
	/*
	 * ==============================================
	 * STEP 1: Initialize Parallel Communication:
	 * 		- Development Status : Done
	 *		- Last modified on: July 23, 2014
	 *						by: Minkwan Kim
	 * =============================================
	 * */

	time_running.initStartTime();

#ifdef	MPI
	MPI_Init(&argc, &argv);					// INITIALIZE MPI
	MPI_Comm_size(MPI_COMM_WORLD, &NP); 	// FINDOUT HOW MANY PROCESSORS IN THERE
	MPI_Comm_rank(MPI_COMM_WORLD, &P);  	// FINDOUT WHICH PROCESSOR I AM
	t0 = MPI_Wtime();
#endif

	if (NP > OP2A_MAX_N_TASK)
		throw Common::ExceptionNPExceed (FromHere(), "Number of processors exceeds MAX_PROCESSOR. Need to adjust value of MAX_N_TASKS");

	kmp_set_defaults("KMP_AFFINITY = scatter");


	/*
	 * =================================================
	 * STEP 2: Show Version information:
	 *		- Development Status : Done
	 *		- Last modified on: July 23, 2014
	 *		-				by: Minkwan Kim
	 * =================================================
	 */
	Version versionOP2A(OP2A_VERSION_MAIN, OP2A_VERSION_SUB, m_now->tm_year + 1900, m_now->tm_mon + 1, m_now->tm_mday, modulename.c_str());

	if (P == 0)
	{
		versionOP2A.info();			// Show the version information
		cout << " --> t = " << time_running.getDeltaT() << "[sec]" << endl;
	}
}
