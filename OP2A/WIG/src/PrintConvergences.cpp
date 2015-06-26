/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 24, 2015
 *      			Author: Minkwan Kim
 *
 * PrintConvergences.cpp
 * 			-  
 *  
 */




#include "Common/include/CodeLocation.hpp"
#include "Common/include/Exception_NPExceed.hpp"
#include "Common/include/StringOps.hpp"

#include "../include/OP2A_Application.hpp"


void ApplicationOP2A::PrintConvergences(bool firststart)
{
	string filename	= problem_setup.name;
	filename += "_convergence.plt";

	ofstream convergence_file;

	if (firststart == true)
	{
		convergence_file.open(filename.c_str());

		// FILE HEADER //
		convergence_file << "TITLE = \""	<< problem_setup.name<< " SOLUTION CONVERGENCES \"" << endl;
		convergence_file << "VARIABLES = \"Iteration number\", \"Residue<sub>max</sub>\", \"Residue<sub>L2</sub>\", \"dt\", \"CFL number\" " << endl;

		//convergence_file << iter << " " << RHS_max << " " << RHS_2 << " " << dt << " " << CFLNumber << endl;
	}
	else
	{
		convergence_file.open(filename.c_str(), std::ofstream::out | std::ofstream::app);
		convergence_file << iter << " " << RHS_max << " " << RHS_2 << " " << dt << " " << CFLNumber << endl;
	}
}
