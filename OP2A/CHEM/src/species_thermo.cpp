/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 15, 2015
 *      			Author: Minkwan Kim
 *
 * species_thermo.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"


// Constructor and Destructor
Species_thermo::Species_thermo()
{
	// THERMODYNAMIC PROPERTIES
	Cv[0]			= 0.0;								/* Translational/Rotational heat capacity */
	Cv[1]			= 0.0;
	R				= 0.0;								/* Gas constant */

	n_vib_lvl		= 0;								/* Number of vibrating mode*/
	n_elec_lvl		= 0;								/* Number of electronic level */

	int i_index;

	for (i_index = 0; i_index <= MAX_ELETRIC_LEVEL-1; i_index++)
	{
		g[i_index]				= 0;					/* Degeneracies of energy levels */
		theta_el[i_index]		= 0.0;			/* Characteristic electronic temperature */
	}

	Z_inf			= 0.0;
	T_star			= 0.0;
	T_ref			= 0.0;
	d_ref			= 0.0;
	omega			= 0.0;

	for (i_index = 0; i_index <= MAX_VIB_MODE-1; i_index++)
	{
		theta_vib[i_index]	= 0.0;				/* Characteristic temp of Vibrational energy */
	}

	int j_index;
	for (i_index = 0; i_index <= 4; i_index++)
		for (j_index = 0; j_index <= 12; j_index++)
			lerc[i_index][j_index]	= 0.0;							/* LeRC curve fits for specific heat, enthalpy and entropy */


	n_vib_ele	= 0;
	for (int i = 0; i <= MAX_V_VIB_RELAX-1; i++)
	{
		kev[i][0]	= 0.0;
		kev[i][1]	= 0.0;
		kev[i][2]	= 0.0;
		kev[i][3]	= 0.0;
	}
}

Species_thermo::~Species_thermo()
{
}




/*
 * =====================================================
 * 	Internal Functions
 * ======================================================
 */

// F01: Assign DATA
void Species_thermo::assign_data(Species_basic species)
{
	// Assign Basic thermodynamic data
	R		= Ru/species.M;		// Gas Constant
	Cv[0]	= 1.5*R;			// Translational specific heat

	switch(species.type)
	{
	case MOLECULE:
		Cv[1]	= R;
		break;
	case ATOM:
		Cv[1]	= 0.0;
		break;
	case ELECTRON:
		Cv[1]	= 0.0;
		break;
	}


	// ROTATIONAL Non-eq DATA
	int 		num;
	double 		Z_inf_temp, T_star_temp, T_ref_temp, omega_temp, d_ref_temp;

	string		name_temp;
	ifstream	database;

	database.open(SPECIES_DATA_ROT_FILE);

	if (species.type == MOLECULE)
	{
		if (database.is_open())
		{
			while (! database.eof())
			{
				database >> name_temp	>> Z_inf_temp >> T_star_temp >> T_ref_temp >> omega_temp >> d_ref_temp;

				num	= max(species.name.size(), name_temp.size());
				if (name_temp.compare(0, num, species.name)	== 0)
				{
					Z_inf	= Z_inf_temp;
					T_star	= T_star_temp;
					T_ref	= T_ref_temp;
					omega	= omega_temp;
					d_ref	= d_ref_temp;
				}
			}
			database.close();
		}
		else
		{
			Error_message_type	error;
			error.location_primary_name 	= "species_thermo.cpp";
			error.location_secondary_name 	= "N/A";
			error.message					= "CANNOT FIND A SPECIES DATABASE FILE FILE. PLEASE CHECK THE FILE";
			error.print_message();
		}
	}



	// VIBRATIONAL Non-eq DATA
	double	theta_vib_temp;

	n_vib_lvl	= 0;
	database.close();
	database.open(SPECIES_DATA_VIB_FILE);

	if (species.type == MOLECULE)
	{
		if (database.is_open())
		{
			while (! database.eof())
			{
				database >> name_temp	>> theta_vib_temp;

				num	= max(species.name.size(), name_temp.size());
				if (name_temp.compare(0, num, species.name)	== 0)
				{
					theta_vib[n_vib_lvl]	= theta_vib_temp;
					n_vib_lvl++;
				}
			}
			database.close();
		}
		else
		{
			Error_message_type	error;
			error.location_primary_name 	= "species_thermo.cpp";
			error.location_secondary_name 	= "N/A";
			error.message					= "CANNOT FIND A SPECIES DATABASE FILE FILE. PLEASE CHECK THE FILE";
			error.print_message();
		}
	}



	// Electronic Energy modes
	int		g_i_temp;
	double	theta_el_temp;

	n_elec_lvl	= 0;
	database.close();
	database.open(SPECIES_DATA_ELE_FILE);


	if (database.is_open())
	{
		while (! database.eof())
		{
			database >> name_temp	>> theta_el_temp	>> g_i_temp;
			num	= max(species.name.size(), name_temp.size());

			if (name_temp.compare(0, num, species.name)	== 0)
			{
				theta_el[n_elec_lvl]	= theta_el_temp;
				g[n_elec_lvl]			= g_i_temp;
				n_elec_lvl++;
			}
		}
		database.close();
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "species_thermo.cpp";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A SPECIES ELECTRONIC STATE DATABASE FILE. PLEASE CHECK THE FILE";
		error.print_message();
	}



	// Enthalpy
	int 		i_index;
	int			num_lerc;
	double		data_temp[13];

	num_lerc	= 0;
	database.open(SPECIES_DATA_LERC_FILE);

	if (database.is_open())
	{
		while (! database.eof())
		{
			database >> name_temp	>> data_temp[0] >> data_temp[1] >> data_temp[2] >> data_temp[3] >> data_temp[4] >> data_temp[5] >> data_temp[6] >> data_temp[7] >> data_temp[8] >> data_temp[9] >> data_temp[10] >> data_temp[11] >> data_temp[12];
			num	= max(species.name.size(), name_temp.size());

			if (name_temp.compare(0, num, species.name)	== 0)
			{
				for (i_index = 0; i_index <= 12; i_index++)
					lerc[num_lerc][i_index]	= data_temp[i_index];

				num_lerc++;
			}
		}

		database.close();
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "species_thermo.cpp";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A SPECIES LeRC DATABASE FILE. PLEASE CHECK THE FILE";
		error.print_message();
	}


	// Vibrational-Electron Relaxation
	int		n_kev_temp;
	database.close();
	database.open(SPECIES_DATA_K_EV_FILE);

	if (database.is_open())
	{
		while (! database.eof())
		{
			database >> name_temp	>> n_kev_temp;
			num	= max(species.name.size(), name_temp.size());

			if (name_temp.compare(0, num, species.name)	== 0)
			{
				n_vib_ele	= n_kev_temp;

				for (int k = 0; k <= n_kev_temp-1; k++)
				{
					database >> kev[k][0] >>  kev[k][1] >>  kev[k][2] >>  kev[k][3];
				}
			}
		}
		database.close();
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "species_thermo.cpp";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A SPECIES VIBRATION-ELECTRON RELAXATION DATABASE FILE. PLEASE CHECK THE FILE";
		error.print_message();
	}

}

