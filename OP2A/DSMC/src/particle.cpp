/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: Dec 12, 2014
 *      			Author: Minkwan Kim
 *
 * particle.cpp
 * 			-  
 *  
 */


#include "../include/particle.hpp"


/*
 * Standard particle type for DSMC
 */
particle_basic_DSMC::particle_basic_DSMC()
{
	species	= -1;

	for (int i = 0; i <= 2; i++)
	{
		X[i]	= 0.0;
		V[i]	= 0.0;
	}

	for (int i = 0; i <= 2; i++)	E[i]	= 0.0;
}

particle_basic_DSMC::~particle_basic_DSMC()
{

}



/*
 * Standard particle type for PIC
 */
particle_basic_PIC::particle_basic_PIC()
{
	W_ratio	= 0.0;
}

particle_basic_PIC::~particle_basic_PIC()
{

}


