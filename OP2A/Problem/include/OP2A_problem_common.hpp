/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_problem_common.hpp
 * 			-  
 *  
 */

#ifndef OP2A_PROBLEM_COMMON_HPP_
#define OP2A_PROBLEM_COMMON_HPP_

#include <string>
#include <vector>



using namespace std;


/*
 * ==================================
 * 	COMMON PROBLEM SETUP CLASS
 * ==================================
 */
class PROBLEM_COMMON
{
public:
	string 			name;							/* PROBLEM NAME				*/

	string			mesh_file_name;					/* Grid file name			*/
	unsigned int	mesh_file_type;					/* Grid file type			*/

	string			output_file_name;				/* Output file name			*/
	unsigned int	output_file_type;				/* Output file type			*/
	double 			grid_factor;
	bool			use_multi_file;					/* Multi file				*/
	unsigned int	itv_save;						/* Interval for writing restart file				*/
	unsigned int	itv_result;						/* Interval for writing result file					*/


	unsigned int	n_current;						// Current simulation time step
	unsigned int	n_total;            			// Maximum number of time steps
	unsigned int	n_startup;        				// Number of time steps until stationary
	double			convergence_criterion;			// CONVERGENCE CRETERION


	bool			is_axisymmetric;				// Axisymmetric simulation
	bool			is_viscous;						// viscous flow

	int				DIM;							// DIMENSION of simulation
	int				NE;								// Number of energy equation
	int				NCM;							// Number of cell in mesh

	int 			NP;									/* Number of CPUs	*/
	int 			NT;									// Number of thread
	int 			P;									// Current CPU ID;

	// Constructor and Destructor
	PROBLEM_COMMON(void);
	~PROBLEM_COMMON(void);


	// M::01 - Read program input file
	void read_COMMON(string file_name);
};



#endif /* OP2A_PROBLEM_COMMON_HPP_ */
