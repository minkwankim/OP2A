/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 15, 2015
 *      			Author: Minkwan Kim
 *
 * specied_trans.cpp
 * 			-  
 *  
 */


#include "../include/species.hpp"
#include "../include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"



// Constructor and Destructor
Species_trans::Species_trans()
{
	Blottner[0]		= 0.0;						// Coefficient for Blottner model: As, Bs, Cs
	Blottner[1]		= 0.0;						// Coefficient for Blottner model: As, Bs, Cs
	Blottner[2]		= 0.0;						// Coefficient for Blottner model: As, Bs, Cs
	Sutherland[0]	= 0.0;						// Coefficient for Sutherland model: mu0, T0, S
	Sutherland[1]	= 0.0;						// Coefficient for Sutherland model: mu0, T0, S
	Sutherland[2]	= 0.0;						// Coefficient for Sutherland model: mu0, T0, S
	sigma			= 0.0;						// Collision diameter in Ansstrong
	T_eps			= 0.0;

	Data_Blottner	= false;
	Data_Sutherland	= false;
	Data_Kinetic	= false;
	Data_Omega		= false;
	Data_Gupta		= false;


	for (int i_index = 0; i_index <= MAX_NS_PROBLEM-1; i_index++)
	{
		CCS[i_index][0][0]	= 0.0;
		CCS[i_index][0][1]	= 0.0;
		CCS[i_index][0][2]	= 0.0;
		CCS[i_index][0][3]	= 0.0;

		CCS[i_index][1][0]	= 0.0;
		CCS[i_index][1][1]	= 0.0;
		CCS[i_index][1][2]	= 0.0;
		CCS[i_index][1][3]	= 0.0;

		for (int j = 0; j <= 11; j++)
		{
			OMEGA[i_index][j][0]	= 0.0;
			OMEGA[i_index][j][1]	= 0.0;
			OMEGA[i_index][j][2]	= 0.0;
		}
	}
}

Species_trans::~Species_trans()
{

}


/*
 * =====================================================
 * 	Internal Functions
 * =====================================================
 */
// F01: Assign Data
void Species_trans::assign_data(Species_basic species)
{
	int			num;
	double		temp1, temp2, temp3;
	string		name_temp;
	ifstream	datafile;


	// 1. Blottner viscosity
	datafile.open(SPECIES_DATA_BLOT_FILE);
	if (datafile.is_open())
	{
		while (! datafile.eof())
		{
			datafile >> name_temp	>> temp1 >> temp2 >> temp3;

			num	= max(species.name.size(), name_temp.size());
			if (name_temp.compare(0, num, species.name)	== 0)
			{
				Blottner[0]		= temp1;
				Blottner[1]		= temp2;
				Blottner[2]		= temp3;
				Data_Blottner	= true;
			}
		}
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "species_trans";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A SPECIES DATABASE FILE FOR BLOTTNER MODEL. PLEASE CHECK THE FILE";
		error.print_message();
	}
	datafile.close();


	// 2. Sutherland model
	datafile.open(SPECIES_DATA_SUTH_FILE);
	if (datafile.is_open())
	{
		while (! datafile.eof())
		{
			datafile >> name_temp	>> temp1 >> temp2 >> temp3;

			num	= max(species.name.size(), name_temp.size());
			if (name_temp.compare(0, num, species.name)	== 0)
			{
				Sutherland[0]	= temp1;
				Sutherland[1]	= temp2;
				Sutherland[2]	= temp3;
				Data_Sutherland	= true;
			}
		}
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "species_trans";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A SPECIES DATABASE FILE FOR SUTHERLAND MODEL. PLEASE CHECK THE FILE";
		error.print_message();
	}
	datafile.close();


	// 3. Kinetic Model
	datafile.open(SPECIES_DATA_KINE_FILE);
	if (datafile.is_open())
	{
		while (! datafile.eof())
		{
			datafile >> name_temp	>> temp1 >> temp2;

			num	= max(species.name.size(), name_temp.size());
			if (name_temp.compare(0, num, species.name)	== 0)
			{
				sigma			= temp1;
				T_eps			= temp2;
				Data_Kinetic	= true;
			}
		}
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "species_trans";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A SPECIES DATABASE FILE FOR KINETIC MODEL. PLEASE CHECK THE FILE";
		error.print_message();
	}
	datafile.close();
}
