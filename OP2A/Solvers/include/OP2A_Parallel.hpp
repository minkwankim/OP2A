/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 15, 2015
 *      			Author: Minkwan Kim
 *
 * Preprocessing.hpp
 * 			-  
 *  
 */
#ifndef PREPROCESSING_HPP_
#define PREPROCESSING_HPP_


#include "Common/include/Common.hpp"
#include "Common/include/OP2A.hpp"
#include "Common/include/Time_StopWatch.hpp"



namespace OP2A{

#define MAX_NUM_THREADS_PER_CPU	2

class OP2A_Parallel
{
public:
	const OP2A::OPuint	NP;		// Total number of processor
	OP2A::OPuint		P;		// Current processor number
	OP2A::OPuint		NT;		// Number of thread

	OP2A_Parallel(bool is_MPI, bool is_OpenMP, bool ) const;

private:
	const bool m_is_MPI;
	const bool m_is_OpenMP;
	const bool m_is_MIC;
};

}


#endif /* PREPROCESSING_HPP_ */
