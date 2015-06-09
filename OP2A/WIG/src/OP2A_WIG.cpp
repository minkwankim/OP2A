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

#include "Common/include/Map1D.hpp"
#include "Common/include/Map2D.hpp"
#include "Common/include/Vector1D.hpp"
#include "Common/include/Vector2D.hpp"

#include "DATA/include/DataStorage.hpp"
#include "DATA/include/DataStorageVector.hpp"

#include "Math/include/OP2A_Vector.hpp"


using namespace OP2A::Setup;
using namespace OP2A;


int main(int argc, char *argv[]) {
	ApplicationOP2A application(OP2A_OPENMP, OP2A_CPU, 23, "OP2A_setup.prob");

	application.preparation(argc, argv, "CFD");


	/*
	 * =======================================================
	 * STEP 3: Read Problem information
	 * 		- Development Status: Version 1.0
	 * 		- Last modified on: July 23, 2014
	 * 						by: Minkwan Kim
	 * @todo need to update using Object Oriented programming method
	 * =======================================================
	 */
	 application.problem_setup.read("Problem_setup_v2.prob");




	 /*
	  * TEST
	  */
	 Map1D <string, int> temp_map(10);
	 temp_map.insert("rho1", 1);
	 temp_map.insert("rho2", 2);
	 temp_map.insert("rho3", 3);
	 temp_map.insert("rho4", 4);
	 Data::DataStorage	data_temp1("V", 4, temp_map);
	 Data::DataStorage	data_temp2("Q", 4, temp_map);

	 Data::DataStorageVector	data_temp(2);
	 data_temp.data[0]	= data_temp1;
	 data_temp.data[1]	= data_temp2;


	 GRID::Grid	grid_test;
	 grid_test.readMeshData("Stardust.msh", GRID::GridDataType::FLUENT);
	 grid_test.processingGridData(1000, false);





	 Math::VECTOR	vector_test1(1.0, 2.0, 3.0);
	 double test = vector_test1(1);

	 int a_test	= vector_test1.dimension();

	 a_test = 1;







/*

	 * =========================================================
	 * STEP 4: Read Species/chemistry data (For NOEQ-CFD mode)
	 * 		- Development Status: Version 1.2a (Need to improve)
	 * 		- Last modified on: July 23, 2014
	 * 						by: Minkwan Kim
	 * =========================================================

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





	 * ======================================================================
	 * STEP 5: GRID GENERATION and/or READ (Unstructured Catersian grid)
	 * 		- Development Status: Version 1.0a
	 * 		- Last modified on: July 23, 2014
	 * 						by: Minkwan Kim
	 * ======================================================================

	 GRID_CLASS 	grid;
	 preprocessing_grid_read_process_ver1(problem.mesh_file_name, problem.mesh_file_type, problem.is_axisymmetric, problem.grid_factor, grid);
	 line_finder(&grid, grid.grid_line.lines, grid.grid_line.lines_bd, grid.grid_line.cell_line_info, grid.grid_line.num_lines);

	 vector < vector<int> >		node_shared_cell_list		= vector_2D<int>(grid.NNM+1, 9, 0);
	 vector < vector<double> >	node_shared_cell_weighting	= vector_2D<double>(grid.NNM+1, 9, 0.0);
	 calculate_node_value_find_make_lists(grid, node_shared_cell_list, node_shared_cell_weighting, problem.NT);





	 * =========================================================
	 * STEP 6: CFD Variable setup
	 * 		- Development Status: Version 1.1a (Need to improve)
	 * 		- Last modified on: July 23, 2014
	 * 						by: Minkwan Kim
	 * =========================================================

	 vector<SOL_CFD>	Solution_data(problem.multi_fluid);
	 switch(problem.multi_fluid)
	 {
	 case 1:
	 Solution_data[0].setup.NS	= species[0].NS;
	 Solution_data[0].setup.ND	= problem_fluid[0].DIM;
	 Solution_data[0].setup.NE	= problem_fluid[0].NE;
	 Solution_data[0].setup.assign_variables(problem_fluid[0].NER, problem_fluid[0].NEV, problem_fluid[0].NEE, fluid_info[0]);
	 break;

	 case 2:
	 Solution_data[0].setup.NS	= species[0].NS;
	 Solution_data[0].setup.ND	= problem_fluid[0].DIM;
	 Solution_data[0].setup.NE	= problem_fluid[0].NE;
	 Solution_data[0].setup.assign_variables(problem_fluid[0].NER, problem_fluid[0].NEV, 0, fluid_info[0]);


	 Solution_data[1].setup.NS	= 1;
	 Solution_data[1].setup.ND	= problem_fluid[0].DIM;
	 Solution_data[1].setup.NE	= 1;
	 Solution_data[1].setup.assign_variables(problem_fluid[0].NER, 0, 0, fluid_info[1]);
	 break;
	 }

	 for (int f = 0; f <= problem.multi_fluid-1; f++)
	 Solution_data[f].allocate_size(grid.NNM, grid.NFM, grid.NCM, grid.NGM, problem_fluid[f].is_viscous, problem_fluid[f].TIME_INTEGRATION_METHOD, problem_fluid[f].is_axisymmetric);






	 * =========================================================
	 * STEP 7: ICs
	 * 		- Development Status: Version 1.1a (Need to improve)
	 * 		- Last modified on: July 23, 2014
	 * 						by: Minkwan Kim
	 * =========================================================

	 vector< vector < vector < double > > > 	rho_IC(problem.multi_fluid);
	 vector< vector < vector < double > > > 	u_IC(problem.multi_fluid);
	 vector < vector < double > > 			T_IC(problem.multi_fluid);
	 vector < vector < double > > 			Tr_IC(problem.multi_fluid);
	 vector < vector < double > > 			Tv_IC(problem.multi_fluid);
	 vector < vector < double > > 			Te_IC(problem.multi_fluid);
	 for (int f = 0; f <= problem.multi_fluid-1; f++)
	 {
	 rho_IC[f]	= problem_fluid[f].IC.rho_s;
	 u_IC[f]		= problem_fluid[f].IC.v;

	 T_IC[f]		= problem_fluid[f].IC.T;
	 Tr_IC[f]	= problem_fluid[f].IC.Tr;
	 Tv_IC[f]	= problem_fluid[f].IC.Tv;
	 Tr_IC[f]	= problem_fluid[f].IC.Tr;
	 }





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



	 * ===================================================================
	 * STEP 10: Preparing Output data
	 * 		- Development Status: Version 1.0
	 * 		- Last modified on: Jan 22, 2015
	 * 						by: Minkwan Kim
	 *====================================================================


	 // Step 8.1 Assign Variable names
	 string	density_pre		= "Density, <greek>r</greek><sub>";
	 string	density_end		= "</sub> [kg/m<sup>-3</sup>]";
	 string	velocity_pre	= "veclocity, ";
	 string	velocity_end	= " [m/s]";
	 string	temperature_pre	= "Temperature, T<sub>";
	 string	temperature_end	= "</sub> [K]";

	 vector < vector <string> > variable_names(problem.multi_fluid);
	 for (int f = 0; f <= problem.multi_fluid-1; f++)
	 {
	 variable_names[f].resize(Solution_data[f].setup.VAR);

	 for (int s = 0; s <= Solution_data[f].setup.NS-1; s++)
	 {
	 int s_ID			= species[f].whereis[s];
	 variable_names[f][s]	= density_pre + species_entire[s_ID].basic_data.name + density_end;
	 }

	 variable_names[f][Solution_data[f].setup.NS]	= velocity_pre + "u" + velocity_end;
	 variable_names[f][Solution_data[f].setup.NS+1]	= velocity_pre + "v" + velocity_end;
	 if (Solution_data[f].setup.ND == 3)	variable_names[f][Solution_data[f].setup.NS+2]	= velocity_pre + "w" + velocity_end;

	 variable_names[f][Solution_data[f].setup.NS + Solution_data[f].setup.ND]	= temperature_pre + "tra" + temperature_end;
	 if (Solution_data[f].setup.ID_T[ROT] != 0)	variable_names[f][Solution_data[f].setup.NS+Solution_data[f].setup.ND+Solution_data[f].setup.ID_T[ROT]]	= temperature_pre + "rot" 	+ temperature_end;
	 if (Solution_data[f].setup.ID_T[VIB] != 0)	variable_names[f][Solution_data[f].setup.NS+Solution_data[f].setup.ND+Solution_data[f].setup.ID_T[VIB]]	= temperature_pre + "vib" 	+ temperature_end;
	 if (Solution_data[f].setup.ID_T[ELE] != 0)	variable_names[f][Solution_data[f].setup.NS+Solution_data[f].setup.ND+Solution_data[f].setup.ID_T[ELE]]	= temperature_pre + "e" 	+ temperature_end;
	 }

	 #ifdef MPI
	 if (P == 0)	cout << "  Write initial Output file  [t = " << MPI_Wtime()-t0 << endl;
	 #else
	 cout << "  Write initial Output file" << endl;
	 #endif

	 // Step 8.2 Write Initial data
	 int	time_strand_ID = 0;
	 vector < vector < vector <double> > > V_print(problem.multi_fluid);
	 for (int f = 0; f <= problem.multi_fluid-1; f++)	V_print[f]	= Solution_data[f].Vc;
	 OP2A_data_print_tecplot_multi_ver1(P, grid, V_print, variable_names, problem.name, problem.output_file_name, problem.multi_fluid);
	 OP2A_data_print_restart(P, grid, V_print, "restart.dat", problem.multi_fluid, problem.n_current, 0);

	 #ifdef MPI
	 if (P == 0)	cout << "  --->Done  [t = " << MPI_Wtime()-t0 << endl;
	 #else
	 cout << "  --->Done" << endl;
	 #endif



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

