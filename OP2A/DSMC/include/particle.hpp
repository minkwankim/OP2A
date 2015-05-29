/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: Dec 12, 2014
 *      			Author: Minkwan Kim
 *
 * particle.hpp
 * 			-  
 *  
 */

#ifndef PARTICLE_HPP_
#define PARTICLE_HPP_



class particle_basic_DSMC
{
public:
	int	species;		// Species info

	double X[3];		// Positions
	double V[3];		// Velocities
	double E[3];		// Internal energies


	/*
	 * Constructor and destructors
	 */
	particle_basic_DSMC();
	~particle_basic_DSMC();
};

class particle_basic_PIC : public particle_basic_DSMC
{
public:
	double W_ratio;		// Weight ratio

	/*
	 * Constructor and destructors
	 */
	particle_basic_PIC();
	~particle_basic_PIC();
};

#define particle_type	particle_basic_DSMC

#ifdef PIC
#undef	particle_type
#define particle_type	particle_basic_PIC
#endif


#endif /* PARTICLE_HPP_ */
