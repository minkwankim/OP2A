/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 15, 2015
 *      			Author: Minkwan Kim
 *
 * species_basic.cpp
 * 			-  
 *  
 */

#include "../include/species.hpp"
#include "../include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"

// Constructor and Destructor
Species_basic::Species_basic()
{
	ID				= -1;
	type			= 0;							/* SPECIES TYPE: */
	M				= 0.0;							/* ATOMIC MASS, AMU */
	m				= 0.0;							/* ATOMIC MASS, kg */
	r				= 0.0; 							/* SPECIES RADIUS */
	h0				= 0.0;							/* ENTHALPHY OF FORMATION,[J/kg] */
	Ds				= 0.0;							/* Dissociation energy [J/kg] */
	I				= 0.0;							/* Ionization energy [J/kg] */
	q				= 0.0;							/* SPECIES CHARGE */
	recombine_to	= 0;							/* Recombination species */
}

Species_basic::~Species_basic()
{

}



/*
 * ===============================================================================
 * 		Internal Functions
 * ===============================================================================
 */

// F01: Read data
void Species_basic::read()
{
	// Basic information
	string 		name_temp;
	int			type_temp;								/* SPECIES TYPE: */
	int 		q_temp;									/* SPECIES CHARGE */
	double		M_temp;									/* ATOMIC MASS, AMU */
	double		r_temp; 								/* SPECIES RADIUS */
	double		h0_temp;								/* ENTHALPHY OF FORMATION,[J/kg] */
	double		Ds_temp;								/* Dissociation energy [J/kg] */
	double		I_temp;									/* Ionization energy [J/kg] */
	string		recombination_temp;

	int num;
	string		line;
	ifstream	database;

	database.open(SPECIES_DATA_BASE_FILE);
	if (database.is_open())
	{
		while (! database.eof())
		{
			//getline(database, line);
			//read_line(line);
			database >> name_temp	>> M_temp	>> h0_temp	>> Ds_temp	>> I_temp	>> q_temp	>> type_temp	>> recombination_temp;

			num	= max(name.size(), name_temp.size());
			if (name_temp.compare(0, num, name)	== 0)
			{
				M				= M_temp;
				m				= M * AMU_SI;
				h0				= h0_temp;
				Ds				= Ds_temp;
				I				= I_temp;
				q				= q_temp;
				type			= type_temp;
				recombination	= recombination_temp;
			}
		}

		database.close();
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "species_trans";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A SPECIES DATABASE FILE. PLEASE CHECK THE FILE";
		error.print_message();
	}
}


