/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_WIG.cpp
 * 			-  
 *  
 */

/*
 * ===============================
 * 		Include Header files
 * ===============================
 */
#include <iostream>
#include <fstream>
#include <time.h>
#include <omp.h>
#include <iomanip>


#include "../include/OP2A_Application.hpp"
#include "../include/OP2A_Problem.hpp"

#include "Setup/include/SetupFileReader.hpp"
#include "GRID/include/Grid.hpp"
#include "GRID/include/PrintResult.hpp"


#include "Common/include/Map1D.hpp"
#include "Common/include/Map2D.hpp"
#include "Common/include/Vector1D.hpp"
#include "Common/include/Vector2D.hpp"

#include "CFD/include/VariableConstants.hpp"

#include "DATA/include/DataStorage.hpp"
#include "DATA/include/DataStorage2D.hpp"
#include "DATA/include/DataStorageVector.hpp"

#include "CHEM/include/SpeciesSet.hpp"
#include "CHEM/include/Reaction.hpp"

#include "Math/include/OP2A_Vector.hpp"
#include "Math/include/OP2A_Matrix.hpp"


using namespace OP2A::Setup;
using namespace OP2A;




int main(int argc, char *argv[])
{
	ApplicationOP2A application(OP2A_OPENMP, OP2A_CPU, 23, "OP2A_setup.prob");
	application.preparation(argc, argv, "CFD");


	/*
	 * =======================================================
	 * STEP 1: Read Problem information
	 * 		- Development Status: Version 1.0
	 * 		- Last modified on: July 23, 2014
	 * 						by: Minkwan Kim
	 * @todo need to update using Object Oriented programming method
	 * =======================================================
	 */
	application.show_starting_task("Read Problem Information");
	application.problem_setup.read("Problem_setup_v2.prob");
	application.check_elapsed_time("Read Problem Information");




	/*
	 * =========================================================
	 * STEP 2: Read Species/chemistry data (For NOEQ-CFD mode)
	 * 		- Development Status: Version 1.2a (Need to improve)
	 * 		- Last modified on: July 23, 2014
	 * 						by: Minkwan Kim
	 * =========================================================
	 */
	application.show_starting_task("Read Species/Chemisty Data");
	application.preprocessing_species();
	/*
	 * @todo need to add Chemical reaction data reader
	 */
	CHEM::Reaction testreaction("N+N2=2N+N+N");
	testreaction.data_completing(application.species_set.species);







	/* ======================================================================
	 * STEP 3: GRID GENERATION and/or READ (Unstructured Catersian grid)
	 * 		- Development Status: Version 1.0a
	 * 		- Last modified on: June 12, 2015
	 * 						by: Minkwan Kim
	 * ======================================================================
	 */
	application.show_starting_task("Read/Generate Grid and allocate solution Data");
	application.preprocessing_grid();




	/*
	 * =======================================================================
	 * STEP 4: INITIALIZE FLOW CONDITION/DATA
	 * 		- Development Status: Version 1.0
	 * 		- Last modified on: June 22, 2015
	 * 						by: Minkwan Kim
	 * ======================================================================
	 */
	application.show_starting_task("Initializing Flow Data");
	application.InitializeData(application.problem_setup.IC.INITIALIZE_METHOD, true);




	/*
	 * =====================================================================
	 * STEP 5: Apply Boundary conditions
	 * 		- Development Status: Version 1.0
	 * 		- Last modified on: Jan 22, 2015
	 * 						by: Minkwan Kim
	 * ======================================================================
	 */
	application.show_starting_task("Applying Inviscid Boundary Condition");
	application.ApplyBCInviscidNormal();
	application.check_elapsed_time("Applying Inviscid Boundary Condition");

	if (application.problem_setup.is_viscous == true)
	{
		application.show_starting_task("Applying Viscous Boundary Condition");
		application.ApplyBCViscousNormal();
		application.check_elapsed_time("Applying Viscous Boundary Condition");
	}
	/*
	 * @todo Need to add transport properties calculator
	 */


	/*
	 * ===================================================================
	 * STEP 6: Preparing Output data
	 * 		- Development Status: Version 1.0
	 * 		- Last modified on: June 15, 2015
	 * 						by: Minkwan Kim
	 *====================================================================
	 */
	// Result Data
	application.show_starting_task("Print solution Data");
	application.print_result(string(NAME_V), 1);

	// Restart Data
	application.show_starting_task("Save restart Data");
	application.print_restartCFD(string(NAME_V));








	/*
	 * =====================================================================
	 * Physical Solving Loop: CFD Part
	 *
	 * @author	Minkwan Kim
	 * @version 1.0 June/24/2015
	 * ======================================================================
	 */
	application.PrintConvergences(true);
	application.show_starting_task("Start OP2A-WIG Module....");

	while (application.problem_setup.n_total && application.termination != true)
	{
		/*
		 * 1. Calculate time step and CFL number
		 */
		application.CalcualteCFL();
		application.Calcualtedt();
		application.t_simulation	+= application.dt;


		/*
		 * 2. Inviscid Part
		 */
		application.ApplyBCInviscidNormal();
		application.CalculateFluxInviscid();
		application.CalculateResidueInviscid();


		/*
		 * 3. Calculate Source Terms
		 */
		application.CalculateSourceTerm();



		/*
		 * 6. Calculate Residue Norms
		 */
		application.CalcualtedResidueNorms();


		/*
		 * 7. Time integral and Update
		 */
		application.TimeIntegrate();
		if(application.P == 0)
		{
			cout << "[Iteration]= " << application.iter << scientific << setprecision(8)
									<< "   [Max. Residual]= " 	<< application.RHS_max
									<< "   [L2 Residual]= " 	<< application.RHS_2
									<< "   [dt]= " << scientific << application.dt
									<< "   [cfl]=" << application.CFLNumber << endl;
		}
		application.PrintConvergences(false);
		application.iter++;

		/*
		 * 8. Print result
		 */
		// Result Data
		if (application.iter % application.problem_setup.itv_result == 0)
		{
			application.show_starting_task("Print solution Data");
			application.print_result(string(NAME_V), 1);

			application.show_starting_task("Save restart Data");
			application.print_restartCFD(string(NAME_V));
		}

		// Restart Data
		if (application.iter % application.problem_setup.itv_save == 0)
		{
			application.show_starting_task("Save restart Data");
			application.print_restartCFD(string(NAME_V));
		}


	}

	cout << "HAHAHHA END" << endl;
	return (0);
}

