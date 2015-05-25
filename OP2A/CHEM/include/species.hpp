/*
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: May 27, 2014
 *      			Author: Minkwan Kim
 *
 * species.hpp
 * 			-  
 *  
 */

#ifndef SPECIES_HPP_
#define SPECIES_HPP_


#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include "OP2A_CHEM_constant.hpp"
#include "Problem/include/OP2A_setup_constants.hpp"

using namespace std;


/*	================================================================
 *		Class NUM 1
 * 			- Species Basic data
 *  ================================================================
 */
class Species_basic
{
public:
	// SPECIES INDICATOR IN CODE
	int		ID;
	string 	name;

	// Basic information
	int		type;								/* SPECIES TYPE: */
	double 	q;									/* SPECIES CHARGE */
	int 	recombine_to;						/* Recombination species */
	string	recombination;

	double	M;									/* ATOMIC MASS, AMU */
	double	m;									/* ATOMIC MASS, kg */
	double	r; 									/* SPECIES RADIUS */

	double	h0;									/* ENTHALPHY OF FORMATION,[J/kg] */
	double	Ds;									/* Dissociation energy [J/kg] */
	double	I;									/* Ionization energy [J/kg] */

	Species_basic();
	~Species_basic();

	void read();
};





/*
 * ==========================================================================
 * 		Class NUM 2
 * 			- Species Thermodynamic data
 * ==========================================================================
 */
class Species_thermo
{
public:
	// THERMODYNAMIC PROPERTIES
	double	Cv[2];								/* Translational/Rotational heat capacity */
	double	R;									/* Species gas constant */


	// Rotational Energy
	double	Z_inf;
	double	T_star;
	double	T_ref;
	double	omega;
	double	d_ref;


	// Vibrational Energy
	int 	n_vib_lvl;							/* Number of vibrating mode*/
	double	theta_vib[MAX_VIB_MODE];			/* Characteristic temp of Vibrational energy */


	// Electronic Energy
	int 	n_elec_lvl;							/* Number of electronic level */
	int 	g[MAX_ELETRIC_LEVEL];				/* Degeneracies of energy levels */
	double	theta_el[MAX_ELETRIC_LEVEL];		/* Characteristic electronic temperature */


	// Enthalpy
	double	lerc[5][13];							/* LeRC curve fits for specific heat, enthalpy and entropy */

	// Vib-electron relaxation mode
	int		n_vib_ele;
	double kev[MAX_V_VIB_RELAX][4];

	Species_thermo();
	~Species_thermo();

	void assign_data(Species_basic species);
};





/*
 * =======================================================
 * 	Class NUM 3
 * 		- Species Transport Properties data
 * =======================================================
 */
class Species_trans
{
public:
	// Transport properties
	double	Blottner[3];						// Coefficient for Blottner model: As, Bs, Cs
	double	Sutherland[3];						// Coefficient for Sutherland model: mu0, T0, S
	double	sigma;								// Collision diameter in Ansstrong
	double	T_eps;								// Molecular parameter
	double	CCS[MAX_NS_PROBLEM][2][4];				// Collision cross section data for neutral collision
	double	OMEGA[MAX_NS_PROBLEM][12][3];			// Collision cross section data for neutral collision
	double	OMEGA_temp[MAX_NS_PROBLEM][12][3];		// Collision cross section data for neutral collision

	bool	Data_Blottner;
	bool	Data_Sutherland;
	bool	Data_Kinetic;
	bool	Data_Gupta;
	bool	Data_Omega;

	Species_trans();
	~Species_trans();

	void assign_data(Species_basic species);
};





/*
 * ============================================================
 * Class Num 3a
 * 	- Additional variables for Thermal Noneq
 * ============================================================
 */
class Species_NONEQ
{
public:
	double Cv[2];

	Species_NONEQ();
	~Species_NONEQ();

	void assign_Cv_prime(double Cv_tra, double Cv_rot, int NONEQ_ROT);
};





/*
 * Class NUM 4
 * 	- Species data
 */
class Species
{
public:
	Species_basic 	basic;
	Species_thermo	thermo;
	Species_trans	transport;
	Species_NONEQ	noneq;


	Species(void);
	~Species(void);

	void asg_data(string name_sp);
	double Cv_vib(double T);
	double Cv_ele(double T);
	double Cv_VE(double T, int NONEQ_VIB, int NONEQ_E);

	double e_vib(double T);
	double e_ele(double T);
	double e_VE(double T, int NONEQ_VIB, int NONEQ_E);

	double e_internal(double T, double Tr, double Tv, double Te, int NONEQ_ROT, int NONEQ_VIB, int NONEQ_E);
	double h(double T, double Tr, double Tv, double Te, int NONEQ_ROT, int NONEQ_VIB, int NONEQ_E);

};



#define SPECIES_BASIC	Species_basic
#define SPECIES_THERMO	Species_thermo
#define SPECIES_TRANS	Species_trans



class SPECIES_ver2
{
public:
	SPECIES_BASIC	basic_data;
	SPECIES_THERMO	thermodynamic_properties;
	SPECIES_TRANS	transport_properties;


	SPECIES_ver2(void);
	~SPECIES_ver2(void);

	void asg_data(string name_sp);

	double Cv_vib(double T);
	double Cv_ele(double T);
	double Cv_VE(double T);
	double Cv_VEE(double T);

	double e_vib(double T);
	double e_ele(double T);
	double e_VE(double T);
	double e_VEE(double T);

	double e_internal(double T, double Tr, double Tv, double Tele);
	double h(double T, double Tr, double Tv, double Tele);

	double tau_e_vib(double Te, double ne);
};







/*
 * Class NUM 5
 * 	- Species data information
 */
class Species_info
{
public:
	int NS;
	int electron_ID;

	// Summary of species Info
	int numMolecules;
	int moleculesID[MAX_NS_PROBLEM];

	Species_info();
	~Species_info();

	void assign_data(Species *species, int NS_setup, int showflag);
	void assign_data(vector<Species> species, int NS_setup, int showflag);

};



#endif /* SPECIES_HPP_ */
