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

#include "Math/include/OP2A_Vector.hpp"


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





	/*
	 * ===================================================================
	 * STEP 10: Preparing Output data
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


	 vector<SPECIES>	species_entire;
	 read_species_data_set(species_entire, problem.species_file, problem.NS);

	 vector<SPECIES_DATA_BASIC>	species(problem.multi_fluid);
	 vector<int> fluid_info(problem.multi_fluid, 0);

	 switch(problem.multi_fluid)
	 {
	 case 1:
	 species[0].allocate_data(species_entire, DATA_ENTIRE);
	 fluid_info[0]	= 0;
	 break;
	 case 2:
	 species[0].allocate_data(species_entire, DATA_HEAVYSPECIES);
	 species[1].allocate_data(species_entire, DATA_ELECTRON);

	 fluid_info[0]	= 0;
	 fluid_info[1]	= -1;
	 break;
	 }

	 // Update Problem setup for each fluids
	 vector<OP2A_PROBLEM> problem_fluid(problem.multi_fluid);
	 for (int f = 0; f <= problem.multi_fluid-1; f++)	problem_fluid[f]	= problem;
	 for (int f = 0; f <=problem.multi_fluid-1; f++)		problem_fluid[f].adjust(species[f]);
	 for (int f = 0; f <=problem.multi_fluid-1; f++)		problem_fluid[f].NS	= species[f].NS;

	 REACTION_DATA_ver2 reactions;
	 reactions.NS			= species_entire.size();
	 reactions.species_data	= species_entire;
	 if (problem.NONEQ_Chem == true)	reactions.read_reaction_data(problem.reaction_file);













	 * =======================================================================
	 * STEP 8: INITIALIZE FLOW CONDITION/DATA
	 * 		- Development Status: Version 1.0
	 * 		- Last modified on: Jan 22, 2015
	 * 						by: Minkwan Kim
	 * ======================================================================

	 #ifdef MPI
	 if (P == 0)	cout << endl << "  Initialize flow conditions [t = " << MPI_Wtime()-t0 << endl;
	 #else
	 cout << endl <<"  Initialize flow conditions..." << endl;
	 #endif
	 vector < vector < vector <double> > > Q_IC(problem.multi_fluid);
	 vector < vector < vector <double> > > V_IC(problem.multi_fluid);
	 Calculate_IC_multi_fluid_ver2(Solution_data, rho_IC, u_IC, T_IC, Tr_IC, Tv_IC, Te_IC, Q_IC, V_IC, species, problem.multi_fluid);

	 // Initialize conditions
	 Initialize_flows_multi_fluids(Solution_data, Q_IC, V_IC, species, grid.NCM, problem.IC.INITIALIZE_METHOD, problem.multi_fluid, problem.NT);
	 for (int f = 0; f <= problem.multi_fluid-1; f++)
	 {
	 #pragma omp parallel for num_threads(NT)
	 for (int c = 0; c <= grid.NCM-1; c++)
	 {
	 Solution_data[f].mixture_data_c[c].calculate_data(problem_fluid[f].NS, problem_fluid[f].NER, Solution_data[f].Qc[c], species[f].Cv_t, species[f].Cv_r, species[f].R, species[f].M);
	 CFD_Q_to_V(Solution_data[f].Qc[c], Solution_data[f].Vc[c], Solution_data[f].setup, Solution_data[f].mixture_data_c[c], species[f]);
	 CFD_V_to_W(Solution_data[f].Vc[c], Solution_data[f].Wc[c], Solution_data[f].setup, Solution_data[f].mixture_data_c[c], species[f]);
	 }
	 }

	 if (problem.is_viscous	== true)
	 {
	 for (int f = 0; f <= problem.multi_fluid-1; f++)
	 {
	 Solution_data[f].setup_transport.Le					=	problem_fluid[f].Le;
	 Solution_data[f].setup_transport.conductivity_model	=	problem_fluid[f].model_conductivity;
	 Solution_data[f].setup_transport.viscosity_model	= 	problem_fluid[f].model_viscosity;
	 Solution_data[f].setup_transport.mixing_rule		=	problem_fluid[f].model_mixing_rule;
	 }

	 for (int f = 0; f <= problem.multi_fluid-1; f++)
	 {
	 for (int c = 0; c <= grid.NCM-1; c++)
	 {
	 Solution_data[f].transport_prop_c[c].calculate_transport_properties(Solution_data[f].setup, Solution_data[f].Vc[c], Solution_data[f].mixture_data_c[c].Xs, species[f], Solution_data[f].setup_transport);
	 }
	 }
	 }

	 #ifdef MPI
	 if (P == 0)	cout << "  --->Done  [t = " << MPI_Wtime()-t0 << endl;
	 #else
	 cout << "  --->Done" << endl;
	 #endif





	 * =====================================================================
	 * STEP 9: Apply Boundary conditions
	 * 		- Development Status: Version 1.0
	 * 		- Last modified on: Jan 22, 2015
	 * 						by: Minkwan Kim
	 * ======================================================================

	 #ifdef MPI
	 if (P == 0)	cout << endl << "  Applying boundary conditions [t = " << MPI_Wtime()-t0 << endl;
	 #else
	 cout << endl <<"  Applying boundary conditions..." << endl;
	 #endif

	 //Inviscid case
	 for (int f = 0; f <= problem.multi_fluid-1; f++)
	 CFD_assign_BC_inviscid_complete(Solution_data[f], grid, Q_IC[f], V_IC[f], species[f], NT);
	 cout <<"   - [Step 1] Inviscid BC --> DONE" << endl;

	 // Viscous CASE
	 if (problem.is_viscous	== true)
	 {
	 for (int f = 0; f <= problem.multi_fluid-1; f++)
	 {

	 CFD_assign_BC_viscous(variable_setup[f], grid,
	 Cell_solution_Q[f], 		Cell_solution_V[f], 		Cell_solution_W[f], 			Cell_CFD_data[f], transport_data_c[f],
	 Cell_ghost_solution_Q[f], 	Cell_ghost_solution_V[f], 	Cell_ghost_solution_W[f], 		Cell_CFD_data[f], transport_data_gc[f],
	 transport_setup[f], species[f],
	 problem_fluid[f].adiabatic, problem_fluid[f].catalytic, problem_fluid[f].radiative, problem_fluid[f].Ys, problem_fluid[f].Tw,
	 problem_fluid[f].use_emissivity, problem_fluid[f].emissivity_Tref, problem_fluid[f].emissivity_below, problem_fluid[f].emissivity_above,
	 problem_fluid[f].NT);

	 }
	 }

	 #ifdef MPI
	 if (P == 0)	cout << "  --->Done  [t = " << MPI_Wtime()-t0 << endl;
	 #else
	 cout << "  --->Done" << endl;
	 #endif
*====================================================================




	 * =======================================================
	 * Physical Solving Loop (CFD)
	 * 		- Development Status: Version 1.0
	 * 		- Last modified on: Dec/15/2014
	 * 						by: Minkwan Kim
	 * =======================================================

	 // [PRESETP]: Setting variables
	 int 						n_RHS;
	 double 						RHS_max, 		RHS2;

	 double 						perf, perf2;
	 vector<double>				Characteristic_length(grid.NCM, 0.0);
	 vector <vector<double> >	Velocity	= vector_2D(problem.multi_fluid, grid.NCM, 0.0);

	 #pragma omp parallel for num_threads(problem.NT)
	 for (int i = 0; i <= grid.NCM-1; i++)	Characteristic_length[i]	= grid.cells.data_ptr[i]->characteristic_length;

	 vector<vector <double *> > 			rho_s_ALL	= NEW_vector_2D_ptr(grid.NCM, reactions.NS);
	 vector<vector<vector <double *> > > Ts_ALL		= NEW_vector_3D_ptr(grid.NCM, reactions.NS, 4);
	 OP2A_CFD_assign_rho_and_T_for_all_fluid(grid.NCM, Solution_data, species, rho_s_ALL, Ts_ALL, problem.multi_fluid);

	 vector<vector <double *> > 			S_chem_ALL	= NEW_vector_2D_ptr(grid.NCM, reactions.NS);
	 for (int f = 0; f <= problem.multi_fluid-1; f++)
	 {
	 for (int s = 0; s <= Solution_data[f].setup.NS-1; s++)
	 {
	 int s_ptr	= species[f].whereis[s];
	 for (int c = 0; c <= grid.NCM-1; c++)	S_chem_ALL[c][s_ptr]	= &Solution_data[f].S_source[c][s];
	 }
	 }

	 vector<vector <double> >	kf	= vector_2D(grid.NCM, reactions.NR, 0.0);
	 vector<vector <double> >	kb	= vector_2D(grid.NCM, reactions.NR, 0.0);
	 vector<vector <double> >	Rf	= vector_2D(grid.NCM, reactions.NR, 0.0);
	 vector<vector <double> >	Rb	= vector_2D(grid.NCM, reactions.NR, 0.0);


	 if (P	== 0)	cout << "Start OP2A-WIG...." << endl;


	 while (problem.n_current < problem.n_total && FLAG_TERMINATION != true)
	 {

	 * ========================================================
	 * 1. Calculate time step
	 * 												- ver 1.0
	 *
	 * 		[Original verion written]
	 * 									- by Minkwan Kim
	 * 									- on 20/Jan/2015
	 * 		[last modified]
	 * 									- by Minkwan Kim
	 * 									- on 27/Jan/2015
	 * ========================================================

	 #pragma omp parallel for num_threads(problem.NT)
	 for (int i = 0; i <= grid.NCM-1; i++)
	 {
	 for (int f = 0; f <= problem.multi_fluid-1; f++)
	 {
	 double U	= 0.0;
	 for (int k = 0; k <= problem_fluid[f].DIM-1; k++)	U += pow(Solution_data[f].Vc[i][problem_fluid[f].NS+k], 2.0);

	 Velocity[f][i]	= sqrt(U);
	 }
	 }

	 double cfl;
	 if (problem.TIME_INTEGRATION_METHOD == 0)
	 {
	 cfl	= problem.CFL_start;
	 }
	 else
	 {
	 if (problem.n_current <= problem.iteration_before_1)
	 cfl =  MK_fn_ver1(0, problem.iteration_before_1, problem.CFL_start, 1, problem.n_current);
	 else
	 cfl =  MK_fn_ver1(problem.iteration_before_1, problem.iteration_before_1*5, 1, problem.CFL_max, problem.n_current);
	 }

	 dt	= calculate_dt_multi(grid.NCM, Characteristic_length, Velocity, problem.max_dt, cfl, problem.multi_fluid);
	 t_simulation	+= dt;



	 * ========================================================
	 * 2. Calculate Required variables (premitive and other variables)
	 * 												- ver 1.0
	 *
	 * 		[Original verion written]
	 * 									- by Minkwan Kim
	 * 									- on 20/Jan/2015
	 * 		[last modified]
	 * 									- by Minkwan Kim
	 * 									- on 13/Feb/2015
	 * ========================================================

	 for (int f = 0; f <= problem.multi_fluid-1; f++)
	 CFD_Calculate_all_required_variables_inviscid(Solution_data[f], grid.NCM, species[f], NT, 0);



	 * ========================================================
	 * 3. Calculate Inviscid Part
	 * 												- ver 1.0
	 *
	 * 		[Original verion written]
	 * 									- by Minkwan Kim
	 * 									- on 20/Jan/2015
	 * 		[last modified]
	 * 									- by Minkwan Kim
	 * 									- on 27/Jan/2015
	 * ========================================================

	 // 3.1. Apply Boundary Conditions
	 for (int f = 0; f <= problem.multi_fluid-1; f++)
	 CFD_assign_BC_inviscid_complete(Solution_data[f], grid, Q_IC[f], V_IC[f], species[f], NT);

	 // 3.2. Calculate flux at faces
	 for (int f = 0; f <= problem.multi_fluid-1; f++)
	 {
	 CFD_Flux_inviscid(grid, Solution_data[f], species[f], problem_fluid[f].NUMERICAL_ORDER, problem_fluid[f].LIMITER, problem_fluid[f].is_axisymmetric, 5.0, 1.0e-5, 0.3, problem_fluid[f].TIME_INTEGRATION_METHOD, NT);
	 CFD_residue_inviscid_ver3(grid, Solution_data[f], problem_fluid[f].is_axisymmetric, NT);
	 }

	 if (problem.TIME_INTEGRATION_METHOD != 0)
	 {
	 for (int f = 0; f <= problem.multi_fluid-1; f++)
	 CFD_assign_BC_inviscid_implicit(Solution_data[f].setup, grid, Solution_data[f].Jacobian_inviscid_plus, Solution_data[f].Jacobian_inviscid_minus, NT);
	 }







	 * ========================================================
	 * 5. Calculate Non-equilibrium Part
	 * 												- ver 1.0
	 *
	 * 		[Original verion written]
	 * 									- by Minkwan Kim
	 * 									- on 20/Jan/2015
	 * 		[last modified]
	 * 									- by Minkwan Kim
	 * 									- on 27/Jan/2015
	 * ========================================================

	 // Chemical rections
	 bool flag_source_term = false;
	 if (problem.NONEQ_Chem == true)
	 {
	 if (problem.TIME_INTEGRATION_METHOD == 0)	CFD_source_chemical_NONEQ_ALL_explicit(grid.NCM, rho_s_ALL, Ts_ALL, reactions, S_chem_ALL, NT);
	 else										CFD_source_chemical_NONEQ_ALL_implicit(grid.NCM, rho_s_ALL, Ts_ALL, reactions, kf, Rf, kb, Rb, S_chem_ALL, NT);

	 //for (int f = 0; f <= problem.multi_fluid-1; f++)	CFD_residue_CHEM_NONEQ(grid, Solution_data[f], problem_fluid[f].is_axisymmetric, NT);
	 flag_source_term = true;
	 }

	 if (problem.NEV == true)
	 {
	 CFD_source_NONEQ_vib(grid, Solution_data[0], species[0], NT);
	 flag_source_term = true;

	 }

	 if (flag_source_term == true)
	 {
	 for (int f = 0; f <= problem.multi_fluid-1; f++)
	 CFD_residue_source(grid, Solution_data[f], problem_fluid[f].is_axisymmetric, NT);
	 }












	 * ========================================================
	 * 6. Residue Norm calculation
	 * 												- ver 1.0
	 *
	 * 		[Original verion written]
	 * 									- by Minkwan Kim
	 * 									- on 20/Jan/2015
	 * 		[last modified]
	 * 									- by Minkwan Kim
	 * 									- on 27/Jan/2015
	 * ========================================================

	 CFD_calc_residue_norms(Solution_data[0].setup.NS, Solution_data[0].setup.ND, Solution_data[0].setup.NE, grid, Solution_data[0].Rn, RHS_max, n_RHS, RHS2, NT);

	 for (int f = 1; f <= problem.multi_fluid-1; f++)
	 {
	 int			n_RHS_temp;
	 double 		RHS_max_temp,	RHS2_temp;

	 CFD_calc_residue_norms(Solution_data[f].setup.NS, Solution_data[f].setup.ND, Solution_data[f].setup.NE, grid, Solution_data[f].Rn, RHS_max_temp, n_RHS_temp, RHS2_temp, NT);

	 if (RHS_max_temp > RHS_max)
	 {
	 RHS_max = RHS_max_temp;
	 n_RHS	= n_RHS_temp;
	 RHS2	= RHS2_temp;
	 }
	 }

	 if (RHS_max <= problem.convergence_criterion && problem.n_current >= 100)	FLAG_TERMINATION = true;
	 else																		FLAG_TERMINATION = false;












	 * ========================================================
	 * 7. Calculate Jacobian of Source terms
	 * 												- ver 1.0
	 *
	 * 		[Original verion written]
	 * 									- by Minkwan Kim
	 * 									- on 20/Jan/2015
	 * 		[last modified]
	 * 									- by Minkwan Kim
	 * 									- on 27/Jan/2015
	 * ========================================================


	 if (problem.TIME_INTEGRATION_METHOD != 0)
	 {
	 for (int f = 0; f <=problem.multi_fluid-1; f++)
	 {
	 #pragma omp parallel for num_threads(NT)
	 for (int c = 0; c <= grid.NCM-1; c++)
	 {
	 CFD_Calculate_dT_dQ(Solution_data[f].setup, Solution_data[f].Qc[c], Solution_data[f].Vc[c], Solution_data[f].mixture_data_c[c], Solution_data[f].dT_dQ[c], species[f]);
	 CFD_Calculate_dp_dQ(Solution_data[f].setup, Solution_data[f].Vc[c], Solution_data[f].mixture_data_c[c], Solution_data[f].dT_dQ[c], Solution_data[f].dp_dQ[c],	species[f]);
	 }
	 }


	 for (int f = 0; f <=problem.multi_fluid-1; f++)
	 {
	 CFD_Jacobian_Source_terms(grid, Solution_data[f], species[f], reactions, kf, Rf, kb, Rb, problem_fluid[f].is_axisymmetric, problem_fluid[f].is_viscous, problem_fluid[f].NONEQ_Chem, NT);
	 }
	 }






	 * ========================================================
	 * 7. Time integral and Update
	 * 												- ver 1.0
	 *
	 * 		[Original verion written]
	 * 									- by Minkwan Kim
	 * 									- on 20/Jan/2015
	 * 		[last modified]
	 * 									- by Minkwan Kim
	 * 									- on 27/Jan/2015
	 * ========================================================


	 for (int f = 0; f <= problem.multi_fluid-1; f++)
	 {
	 CFD_time_integration(grid, Solution_data[f], dt, problem.is_axisymmetric, problem.TIME_INTEGRATION_METHOD, problem.NT);
	 CFD_update_Q(grid, Solution_data[f], NT);

	 #pragma omp parallel for num_threads(NT)
	 for (int c = 0; c <= grid.NCM-1; c++)
	 Solution_data[f].Qc[c]	= Solution_data[f].Q_new[c];
	 }






	 * ========================================================
	 * 8. Print Result
	 * 												- ver 1.0
	 *
	 * 		[Original verion written]
	 * 									- by Minkwan Kim
	 * 									- on 20/Jan/2015
	 * 		[last modified]
	 * 									- by Minkwan Kim
	 * 									- on 27/Jan/2015
	 * ========================================================

	 problem.n_current++;
	 #pragma omp parallel for num_threads(problem.multi_fluid)
	 for (int f = 0; f <= problem.multi_fluid-1; f++)	problem_fluid[f].n_current	= problem.n_current;

	 if(P == 0)	cout << "[Iteration]= " << problem.n_current << scientific << setprecision(8) << "   [Max. Residual]= " << RHS_max << "   [L2 Residual]= " << RHS2 <<"   [dt]= " << scientific << dt << "  [cfl]=" << cfl << endl;

	 if (problem.n_current % problem.itv_result == 0)
	 {
	 #ifdef MPI
	 if (P==0)	cout << " Writing output data [t]=" << MPI_Wtime()-t0 << endl;
	 #else
	 cout << " Writing output data" << endl;
	 #endif
	 for (int f = 0; f <= problem.multi_fluid-1; f++)	CFD_Calculate_all_required_variables_inviscid(Solution_data[f], grid.NCM, species[f], NT, 0);
	 for (int f = 0; f <= problem.multi_fluid-1; f++)	V_print[f]	= Solution_data[f].Vc;

	 OP2A_data_print_tecplot_multi_ver1(P, grid, V_print, variable_names, problem.name, problem.output_file_name, problem.multi_fluid);
	 OP2A_data_print_restart(P, grid, V_print, "restart.dat", problem.multi_fluid, problem.n_current, 0);
	 }







	 }






























	 * ========================================================
	 * 4. Calculate Viscous Part
	 * 												- ver 1.0
	 *
	 * 		[Original verion written]
	 * 									- by Minkwan Kim
	 * 									- on 20/Jan/2015
	 * 		[last modified]
	 * 									- by Minkwan Kim
	 * 									- on 27/Jan/2015
	 * ========================================================

	 if (problem.is_viscous	== true)
	 {
	 for (int f = 0; f <= NUM_FLUID-1; f++)
	 {
	 #pragma omp parallel for num_threads(problem.NT)
	 for (int c = 0; c <= grid.NCM-1; c++)
	 {
	 transport_data_c[f][c].calculate_transport_properties(variable_setup[f], Cell_solution_V[f].data_ptr[c]->data, Cell_CFD_data[f][c].Xs, species[f], transport_setup[f]);
	 }
	 }

	 for (int f = 0; f <= NUM_FLUID-1; f++)
	 {
	 CFD_assign_BC_viscous(variable_setup[f], grid,
	 Cell_solution_Q[f], 		Cell_solution_V[f], 		Cell_solution_W[f], 			Cell_CFD_data[f], transport_data_c[f],
	 Cell_ghost_solution_Q[f], 	Cell_ghost_solution_V[f], 	Cell_ghost_solution_W[f], 		Cell_CFD_data[f], transport_data_gc[f],
	 transport_setup[f], species[f],
	 problem_fluid[f].adiabatic, problem_fluid[f].catalytic, problem_fluid[f].radiative, problem_fluid[f].Ys, problem_fluid[f].Tw,
	 problem_fluid[f].use_emissivity, problem_fluid[f].emissivity_Tref, problem_fluid[f].emissivity_below, problem_fluid[f].emissivity_above,
	 problem_fluid[f].NT);
	 }

	 double vis_relaxation;
	 //vis_relaxation	= MK_fn_ver1(problem.viscous_relaxation_na, problem.viscous_relaxation_nb, problem.viscous_relaxation_min, problem.viscous_relaxation_max, problem.n_current);
	 vis_relaxation = 0.1;

	 CFD_Flux_viscous_multi_fluids_ver2(variable_setup, grid, problem.is_axisymmetric,
	 node_shared_cell_list, node_shared_cell_weighting,
	 Cell_solution_Q, 		Cell_solution_V, 		Cell_solution_W,		Cell_CFD_data, 			transport_data_c,
	 Cell_ghost_solution_Q, 	Cell_ghost_solution_V, 	Cell_ghost_solution_W,	Cell_ghost_CFD_data, 	transport_data_gc,
	 species,
	 flux_viscous_face,  residue_cell,
	 Jacobian_v_plus, Jacobian_v_minus,
	 transport_setup,
	 0, 1, 0, problem.TIME_INTEGRATION_METHOD, problem.adiabatic, problem.catalytic, problem.radiative,
	 NUM_FLUID, problem.NT, vis_relaxation);

	 if (problem.TIME_INTEGRATION_METHOD != 0)
	 {
	 for (int f = 0; f <= NUM_FLUID-1; f++)	CFD_Flux_Jacobian_(grid, Jacobian_plus[f], Jacobian_minus[f], Jacobian_v_plus[f], Jacobian_v_minus[f], 0.01, variable_setup[f].VAR, problem.NT);
	 }
	 }








	 * ========================================================
	 * 5. Calculate Non-equilibrium Part
	 * 												- ver 1.0
	 *
	 * 		[Original verion written]
	 * 									- by Minkwan Kim
	 * 									- on 20/Jan/2015
	 * 		[last modified]
	 * 									- by Minkwan Kim
	 * 									- on 27/Jan/2015
	 * ========================================================

	 // Chemical rections
	 if (problem.NONEQ_Chem == true)
	 {
	 CFD_Chemical_NONEQ_single(variable_setup[0], grid, Cell_solution_V[0], species[0], reactions,	S_chem[0], problem.NT);
	 //CFD_Chemical_NONEQ_multi(variable_setup, grid, Cell_solution_V, species, reactions, S_chem, NUM_FLUID, problem.NT);
	 CFD_residue_CHEM_NONEQ_multi_fluids(variable_setup, problem.is_axisymmetric, grid, residue_cell, S_chem, NUM_FLUID, problem.NT);
	 }













	 }


	 DELETE Allocated Pointers
	 delete[]	Cell_solution_Q;
	 delete[]	Cell_solution_V;
	 delete[]	Cell_solution_W;
	 delete[]	Cell_CFD_data;

	 delete[]	Cell_ghost_solution_Q;
	 delete[]	Cell_ghost_solution_V;
	 delete[]	Cell_ghost_solution_W;
	 delete[]	Cell_ghost_CFD_data;

	 delete[]	flux_inviscid_face;
	 delete[]	residue_cell;
	 delete[]	cell_dQn;
	 delete[]	cell_Qnew;
	 delete[]	data_print;
	 delete[]	S_chem;
	 delete[]	Jacobian_plus;
	 delete[]	Jacobian_minus;
	 delete[]	Jacobian_source;

	 delete[]	flux_viscous_face;
	 delete[]	transport_data_c;
	 delete[]	transport_data_gc;

	 DELETE_vector_2D_ptr(rho_s_ALL, grid.NCM, reactions.NS);
	 DELETE_vector_3D_ptr(Ts_ALL, grid.NCM, reactions.NS, 4);*/
	cout << "HAHAHHA END" << endl;
	return (0);
}

