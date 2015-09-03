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
#include "Common/include/MultiDimension.hpp"

#include "CFD/include/VariableConstants.hpp"

#include "DATA/include/DataStorage.hpp"
#include "DATA/include/DataStorage2D.hpp"
#include "DATA/include/DataStorageVector.hpp"

#include "CHEM/include/SpeciesSet.hpp"
#include "CHEM/include/Reaction.hpp"

#include "Math/include/OP2A_Vector.hpp"
#include "Math/include/OP2A_Matrix.hpp"
#include "Math/include/OP2A_Math.hpp"


#include "CFD/include/GradientCalculation.hpp"

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
	application.show_starting_task("Read Species/Chemistry Data");
	application.preprocessing_species();

	vector<double> Ys_temp(12, 0.0);
	Ys_temp[0]	= 0.78;
	Ys_temp[1]	= 0.22;
	double D12 = application.species_set.ThermalConductivity_e_Gupta(373, 1.0e17, 373, 101300, Ys_temp);


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




	//TESTing GRAD
	for (int c = 1; c <= application.grid.NCM; c++)
	{
		double x = application.grid.cells[c].geo.x[0];
		double y = application.grid.cells[c].geo.x[1];
		double ff = x*x*x + x*y + y*y + exp(x+y);

		application.grid.cells[c].data1D.data[0].data[0] = ff;
	}

	for (int c = 1; c <= application.grid.NGM; c++)
	{
		double x = application.grid.cells_ghost[c].geo.x[0];
		double y = application.grid.cells_ghost[c].geo.x[1];
		double ff = x*x*x + x*y + y*y + exp(x+y);

		application.grid.cells_ghost[c].data1D.data[0].data[0] = ff;
	}

	vector<vector <double> > grad_phi	= Common::vector_2D(application.grid.NFM+1, 2, 0.0);
	CFD::GradientCalculation::GradientCalculation_LSQR(application.grid, 0, 0, 1, grad_phi);


	double error_max = 0.0;
	for (int c = 1; c <= application.grid.NCM; c++)
	{
		double x = application.grid.cells[c].geo.x[0];
		double y = application.grid.cells[c].geo.x[1];

		double fx = 3.0*x*x + y + exp(x+y);
		double fy = x + 2.0*y + exp(x+y);

		double fx_cal = grad_phi[c][0];
		double fy_cal = grad_phi[c][1];


		double error_x = Math::fabs<double>(fx - fx_cal) / fx * 100.0;
		double error_y = Math::fabs<double>(fy - fy_cal) / fy * 100.0;

		error_max = Math::fmax<double>(error_x, error_max);
		error_max = Math::fmax<double>(error_y, error_max);
	}



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
		 * @author Minkwan Kim
		 * @version 1.0 June/24/2015
		 * [Done]: It does not need to be improved
		 */
		application.CalcualteCFL();						// Calculate CFL number based on explicit/implicit methods
		application.Calcualtedt();						// Calculate time step
		application.t_simulation	+= application.dt;



		/*
		 * 2. Inviscid Part
		 * @author Minkwan Kim
		 * @version 1.0 Aug/19/2015
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

