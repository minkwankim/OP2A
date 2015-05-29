/************************************************************************
			Open-source Multi-Physics Solver - ver. 0.0		
					- Plasma physics module 
				
Copyright(c) 2013 MINKWAN KIM												
							
Initial Developed Date: Jan 16, 2014
                    by: Minkwan Kim

Last Modified Date: Jan 16, 2014
					by: Minkwan Kim


constants.hpp
	- Definition of physical constants for plasma physics
              
					Copyright 2013 MINKWAN KIM
*************************************************************************/

#ifndef _PLASMA_CONSTANTS_H
#define _PLASMA_CONSTANTS_H


namespace OP2A
{
	namespace PLAMSA
	{
#define	BOHR_RADIUS			5.29177e-11						// Bohr radius
#define EPS0_SI 			8.854E-12 						// Vacuum permittivity in SI
#define EPS0_CG 			5.772156649015329E-1
#define ELECTRON_CHARGE_CGS 4.803E-10
#define ELECTRON_CHARGE_SI	1.60217646E-19

//	#define C_BOLTZMANN_CGS 1.3806E-16
//	#define C_BOLTZMANN_SI 1.3806E-23
//
//#define N_AVOGADRO_SI 6.022E26
//#define N_AVOGADRO_CGS 6.022E23
//
//#define ELECTRON_CHARGE_CGS 4.803E-10
//#define ELECTRON_CHARGE_SI 1.60217646E-19
//#define ELECTRON_CHARGE_CGS_2 (ELECTRON_CHARGE_CGS*ELECTRON_CHARGE_CGS)
//
//#define ANGSTRONS2_TO_M2 1.E-20
//#define M2_TO_ANGSTRONS2 1.E20
//
//#define Ru 8314.0
//#define Ru_cgs 8.314
//
//#define Pa_TO_ATM 9.8692E-6
//#define ONE_BAR_JCM3 0.1
//#define TORR_TO_PA 133.3224
//
//#define Me 		9.10938188E-31			  		/*electron mass in kg */
//#define Me_cgs	9.10938188E-28					// Electron mass in g
//#define Me_MOLE 5.485669768136001e-04 			/* Electron Molecular weight */
//#define Cv_e 1.5*Ru/ME_MOLE
//
//#define AMU_SI 1.660538782E-27
//#define AMU_CGS 1.660538782E-24
//
//#define STEFAN_BOLTZMANN_SI 5.6704E-8			/*Stefan-Boltzmann constant in SI unit*/
//#define K_TO_eV 8.67343e-5
//
//#ifndef PI
//	#define PI 3.14159265
//#endif
//
//#define PLANCK_SI				6.6262e-34
//#define PLANCK_SI_BAR			1.0546e-34
//#define SPEED_OF_LIGHT_SI		2.9979e8
//#define MU0_SI					(4*PI*1.0e-7)			/* Permeability of free space */
//#define BOHR_RADIUS_SI			5.2918e-11				/* Bohr radius */
//#define ATOMIC_CROSS_SECTION_SI	8.7974e-21				// Atomic cross section, pi a0^2
//#define eV_TO_K					11605
//#define E_in_eV					1.6022e-19				// Energy for 1V


	}



}



#endif
