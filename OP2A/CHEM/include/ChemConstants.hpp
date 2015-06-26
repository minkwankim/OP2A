/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 10, 2015
 *      			Author: Minkwan Kim
 *
 * ChemConstants.hpp
 * 			-  
 *  
 */
#ifndef CHEMCONSTANTS_HPP_
#define CHEMCONSTANTS_HPP_


namespace OP2A{
namespace CHEM{

#define Ru 8314.0
#define Ru_cgs 8.314

#define C_BOLTZMANN_CGS 1.3806E-16
#define C_BOLTZMANN_SI 1.3806E-23

#define N_AVOGADRO_SI 6.022E26
#define N_AVOGADRO_CGS 6.022E23

#define Me 		9.10938188E-31			  		/*electron mass in kg */
#define Me_cgs	9.10938188E-28					// Electron mass in g
#define Me_MOLE 5.485669768136001e-04 			/* Electron Molecular weight */

#define Cv_e 1.5*Ru/ME_MOLE

#define AMU_SI 	1.660538782E-27
#define AMU_CGS 1.660538782E-24

#define STEFAN_BOLTZMANN_SI 5.6704E-8			/*Stefan-Boltzmann constant in SI unit*/
#define EPS0_SI 8.854E-12 						/* Vacuum permittivity in SI*/
#define EPS0_CG 5.772156649015329E-1
#define K_TO_eV 8.67343e-5

#define PLANCK_SI				6.6262e-34
#define PLANCK_SI_BAR			1.0546e-34
#define SPEED_OF_LIGHT_SI		2.9979e8
#define MU0_SI					(4*PI*1.0e-7)			/* Permeability of free space */
#define BOHR_RADIUS_SI			5.2918e-11				/* Bohr radius */
#define ATOMIC_CROSS_SECTION_SI	8.7974e-21				// Atomic cross section, pi a0^2
#define eV_TO_K					11605
#define E_in_eV					1.6022e-19				// Energy for 1V



#define OP2A_SPECIES_DATA_BASE_FILE		"species_database.dat"
#define OP2A_SPECIES_DATA_ROT_FILE		"species_noneq_rot.dat"
#define OP2A_SPECIES_DATA_VIB_FILE		"species_noneq_vib.dat"
#define OP2A_SPECIES_DATA_ELE_FILE		"species_electronic.dat"
#define OP2A_SPECIES_DATA_K_EV_FILE		"species_kev_coeff_theory.dat"
#define OP2A_SPECIES_DATA_LERC_FILE		"species_LeRC.dat"
#define OP2A_SPECIES_DATA_BLOT_FILE		"species_Blottner.dat"
#define OP2A_SPECIES_DATA_SUTH_FILE		"species_Sutherland.dat"
#define OP2A_SPECIES_DATA_KINE_FILE		"species_Kinetic.dat"
#define OP2A_SPECIES_DATA_PARK_FILE		"species_Parker.dat"
#define OP2A_SPECIES_DATA_OMEGA11_FILE	"species_Omega11.dat"
#define OP2A_SPECIES_DATA_OMEGA22_FILE	"species_Omega22.dat"

}
}
#endif /* CHEMCONSTANTS_HPP_ */
