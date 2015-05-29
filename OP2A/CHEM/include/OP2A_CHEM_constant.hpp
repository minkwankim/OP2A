/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: May 27, 2014
 *      			Author: Minkwan Kim
 *
 * constant.hpp
 * 			-  
 *  
 */

#ifndef _CHEM_CONSTANT_HPP_
#define _CHEM_CONSTANT_HPP_



#define SPECIES_DATA_BASE_FILE		"species_database.dat"
#define SPECIES_DATA_ROT_FILE		"species_noneq_rot.dat"
#define SPECIES_DATA_VIB_FILE		"species_noneq_vib.dat"
#define SPECIES_DATA_ELE_FILE		"species_electronic.dat"
#define SPECIES_DATA_K_EV_FILE		"species_kev_coeff_theory.dat"
#define SPECIES_DATA_LERC_FILE		"species_LeRC.dat"
#define SPECIES_DATA_BLOT_FILE		"species_Blottner.dat"
#define SPECIES_DATA_SUTH_FILE		"species_Sutherland.dat"
#define SPECIES_DATA_KINE_FILE		"species_Kinetic.dat"
#define SPECIES_DATA_PARK_FILE		"species_Parker.dat"
#define SPECIES_DATA_OMEGA11_FILE	"species_Omega11.dat"
#define SPECIES_DATA_OMEGA22_FILE	"species_Omega22.dat"


// Physical Limits
#define MAX_ELETRIC_LEVEL	16
#define	MAX_V_VIB_RELAX		10
#define MAX_VIB_MODE		2


// Reaction names
#define	DISSOCIATION					0
#define EXCHANGE						1
#define DISSOCIATIVE_RECOMBINATION		2
#define CHARGE_EXCHANGE					3
#define	ELECTRON_IMPACT_DISSOCIATION	4
#define	ELECTRON_IMPACT_IONIZATION		5

#define	PARK85							85
#define	PARK90							90
#define	PARK94							94
#define	LeRC							0




// Text values
#define TRA					0
#define ROT					1
#define VIB					2
#define ELE					3
#define VE					4
#define VEE					5
#define	TR					6

#define	MOLECULE	1
#define	ATOM		0
#define	ELECTRON	-1

#define	ION			1
#define NEUTRAL		0

#define	DATA_ENTIRE				0
#define	DATA_NEUTRAL			1
#define	DATA_ION				2
#define	DATA_ELECTRON			3
#define	DATA_HEAVYSPECIES		4
#define	DATA_MOLECULES			5
#define	DATA_ATOMS				6


#endif /* CONSTANT_HPP_ */
