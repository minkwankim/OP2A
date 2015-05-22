/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 15, 2015
 *      			Author: Minkwan Kim
 *
 * species_functions.cpp
 * 			-  
 *  
 */



#include "../include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"



void read_species_ver2(vector<SPECIES> &species,	string file_name, int NS)
{
	int s;
	ifstream	species_setting;
	string	name;

	species.resize(NS);

	species_setting.open(file_name.c_str());

	// READ PROBLEM SETUP DATA FILE
	if (species_setting.is_open())
	{
		for (s = 0; s <= NS-1; s++)
		{
			species_setting >> name;
			species[s].asg_data(name);
			species[s].basic_data.ID	= s;
		}
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "species_functions.cpp";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A SPECIES SETUP DATA FILE. PLEASE CHECK THE FILE";
		error.print_message();
	}

	species_setting.close();
}



void asg_recombination_ver2(vector<SPECIES>& species, int NS)
{
	int i_index, s_index;
	int num;

	for (i_index = 0; i_index <= NS-1; i_index++)
	{
		for (s_index = 0; s_index <= NS-1; s_index++)
		{
			num	= max(species[i_index].basic_data.recombination.size(), species[s_index].basic_data.name.size());
			if (species[i_index].basic_data.recombination.compare(0, num, species[s_index].basic_data.name)	== 0)
			{
				species[i_index].basic_data.recombine_to	= s_index;
				break;
			}
		}
	}
}


void read_species_data_set(vector<SPECIES>& species, string file_name, int NS)
{
	// Read Basic data + Transport Properties data
	read_species(species, file_name, NS);
	asg_recombination(species, NS);

	// Read Collisional integral data

}


















/*
 * Species Functions
 * 		ver - 1.0
 */
void FSpecies_read_species(Species *species, string file_name, int NS)
{
	int s;
	ifstream	species_setting;
	string	name;

	species_setting.open(file_name.c_str());

	// READ PROBLEM SETUP DATA FILE
	if (species_setting.is_open())
	{
		for (s = 0; s <= NS-1; s++)
		{
			species_setting >> name;
			species[s].asg_data(name);
		}
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "species_functions.cpp";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A SPECIES SETUP DATA FILE. PLEASE CHECK THE FILE";
		error.print_message();
	}

	species_setting.close();
}


void FSpecies_read_species2(vector<Species>& species,	string file_name, int NS)
{
	int s;
	ifstream	species_setting;
	string	name;

	species_setting.open(file_name.c_str());

	// READ PROBLEM SETUP DATA FILE
	if (species_setting.is_open())
	{
		for (s = 0; s <= NS-1; s++)
		{
			species_setting >> name;
			species[s].asg_data(name);
		}
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "species_functions.cpp";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A SPECIES SETUP DATA FILE. PLEASE CHECK THE FILE";
		error.print_message();
	}

	species_setting.close();
}







void FSpecies_asg_recombination(Species *species, int NS)
{
	int i_index, s_index;
	int num;

	for (i_index = 0; i_index <= NS-1; i_index++)
	{
		for (s_index = 0; s_index <= NS-1; s_index++)
		{
			num	= max(species[i_index].basic.recombination.size(), species[s_index].basic.name.size());
			if (species[i_index].basic.recombination.compare(0, num, species[s_index].basic.name)	== 0)
			{
				species[i_index].basic.recombine_to	= s_index;
				break;
			}
		}
	}
}

void FSpecies_asg_recombination2(vector<Species>& species, int NS)
{
	int i_index, s_index;
	int num;

	for (i_index = 0; i_index <= NS-1; i_index++)
	{
		for (s_index = 0; s_index <= NS-1; s_index++)
		{
			num	= max(species[i_index].basic.recombination.size(), species[s_index].basic.name.size());
			if (species[i_index].basic.recombination.compare(0, num, species[s_index].basic.name)	== 0)
			{
				species[i_index].basic.recombine_to	= s_index;
				break;
			}
		}
	}
}


void Read_species_DATA(Species *species, string file_name, int NS)
{
	// Read Basic data + Transport Properties data
	FSpecies_read_species(species, file_name, NS);
	FSpecies_asg_recombination(species, NS);

	// Read Collisional integral data
	//FAsg_Omega_integral(species, NS);

}

void Read_species_DATA2(vector<Species>& species, string file_name, int NS)
{
	// Read Basic data + Transport Properties data
	FSpecies_read_species2(species, file_name, NS);
	FSpecies_asg_recombination2(species, NS);

	// Read Collisional integral data
	FAsg_Omega_integral2(species, NS);
}



void Read_species_DATA_v3(vector<Species>& species, string file_name, int NS)
{
	// Read Basic data + Transport Properties data
	FSpecies_read_species2(species, file_name, NS);

}







void FAsg_Omega_integral(Species* species, int NS)
{

	int 		i, j, k;
	int 		num;
	string		SP1, SP2;
	double		temp[12];
	double		coeff[12];
	ifstream	datafile;



	// Omega 11	//
	datafile.open(SPECIES_DATA_OMEGA11_FILE);

	if (datafile.is_open())
	{
		while (! datafile.eof())
		{
			datafile >> SP1	>> SP2;
			for (i = 0; i <= 11; i++)	datafile >> temp[i];
			for (i = 0; i <= 11; i++)	datafile >> coeff[i];

			for (i = 0;  i <= NS-1; i++)
			{
				num	= max(species[i].basic.name.size(), SP1.size());
				if (SP1.compare(0, num, species[i].basic.name)	== 0)
				{
					for (j =0; j <= NS-1; j++)
					{
						num	= max(species[j].basic.name.size(), SP2.size());
						if (SP2.compare(0, num, species[j].basic.name)	== 0)
						{
							for (k = 0;  k <= 11; k++)
							{
								species[i].transport.OMEGA[j][k][1]			= coeff[k];
								species[i].transport.OMEGA_temp[j][k][1]	= temp[k];

								species[j].transport.OMEGA[i][k][1]		= coeff[k];
								species[j].transport.OMEGA_temp[i][k][1]	= temp[k];
							}

						}
					}

				}
			}
		}
	}
	datafile.close();


	// Omega 22	//
	datafile.open(SPECIES_DATA_OMEGA22_FILE);

	if (datafile.is_open())
	{
		while (! datafile.eof())
		{
			datafile >> SP1	>> SP2;
			for (i = 0; i <= 11; i++)	datafile >> temp[i];
			for (i = 0; i <= 11; i++)	datafile >> coeff[i];

			for (i = 0;  i <= NS-1; i++)
			{
				num	= max(species[i].basic.name.size(), SP1.size());
				if (SP1.compare(0, num, species[i].basic.name)	== 0)
				{
					for (j =0; j <= NS-1; j++)
					{
						num	= max(species[j].basic.name.size(), SP2.size());
						if (SP2.compare(0, num, species[j].basic.name)	== 0)
						{
							for (k = 0;  k <= 11; k++)
							{
								species[i].transport.OMEGA[j][k][2]		= coeff[k];
								species[i].transport.OMEGA_temp[j][k][2]	= temp[k];

								species[j].transport.OMEGA[i][k][2]		= coeff[k];
								species[j].transport.OMEGA_temp[i][k][2]	= temp[k];
							}
						}
					}
				}
			}
		}
	}
	datafile.close();

	// Check Status
	for (i = 0; i <= NS-1; i++)
	{
		species[i].transport.Data_Omega	= true;
		for (j = 0; j <= NS-1; j++)
		{
			if (species[i].transport.OMEGA[j][0][1] == 0.0 || species[i].transport.OMEGA[j][0][2] == 0.0)
			{
				species[i].transport.Data_Omega	= false;
				break;
			}
		}
	}
}



void FAsg_Omega_integral2(vector<Species>& species, int NS)
{

	int 		i, j, k;
	int 		num;
	string		SP1, SP2;
	double		temp[12];
	double		coeff[12];
	ifstream	datafile;



	// Omega 11	//
	datafile.open(SPECIES_DATA_OMEGA11_FILE);

	if (datafile.is_open())
	{
		while (! datafile.eof())
		{
			datafile >> SP1	>> SP2;
			for (i = 0; i <= 11; i++)	datafile >> temp[i];
			for (i = 0; i <= 11; i++)	datafile >> coeff[i];

			for (i = 0;  i <= NS-1; i++)
			{
				num	= max(species[i].basic.name.size(), SP1.size());
				if (SP1.compare(0, num, species[i].basic.name)	== 0)
				{
					for (j =0; j <= NS-1; j++)
					{
						num	= max(species[j].basic.name.size(), SP2.size());
						if (SP2.compare(0, num, species[j].basic.name)	== 0)
						{
							for (k = 0;  k <= 11; k++)
							{
								species[i].transport.OMEGA[j][k][1]			= coeff[k];
								species[i].transport.OMEGA_temp[j][k][1]	= temp[k];

								species[j].transport.OMEGA[i][k][1]		= coeff[k];
								species[j].transport.OMEGA_temp[i][k][1]	= temp[k];
							}

						}
					}

				}
			}
		}
	}
	datafile.close();


	// Omega 22	//
	datafile.open(SPECIES_DATA_OMEGA22_FILE);

	if (datafile.is_open())
	{
		while (! datafile.eof())
		{
			datafile >> SP1	>> SP2;
			for (i = 0; i <= 11; i++)	datafile >> temp[i];
			for (i = 0; i <= 11; i++)	datafile >> coeff[i];

			for (i = 0;  i <= NS-1; i++)
			{
				num	= max(species[i].basic.name.size(), SP1.size());
				if (SP1.compare(0, num, species[i].basic.name)	== 0)
				{
					for (j =0; j <= NS-1; j++)
					{
						num	= max(species[j].basic.name.size(), SP2.size());
						if (SP2.compare(0, num, species[j].basic.name)	== 0)
						{
							for (k = 0;  k <= 11; k++)
							{
								species[i].transport.OMEGA[j][k][2]		= coeff[k];
								species[i].transport.OMEGA_temp[j][k][2]	= temp[k];

								species[j].transport.OMEGA[i][k][2]		= coeff[k];
								species[j].transport.OMEGA_temp[i][k][2]	= temp[k];
							}
						}
					}
				}
			}
		}
	}
	datafile.close();

	// Check Status
	for (i = 0; i <= NS-1; i++)
	{
		species[i].transport.Data_Omega	= true;
		for (j = 0; j <= NS-1; j++)
		{
			if (species[i].transport.OMEGA[j][0][1] == 0.0 || species[i].transport.OMEGA[j][0][2] == 0.0)
			{
				species[i].transport.Data_Omega	= false;
				break;
			}
		}
	}
}






void FSpecies_updata_Noneq_data(Species* species, int NS, int NONEQ_ROT, int NONEQ_VIB, int NONEQ_E)
{
	int	i;

	for (i = 0; i <= NS-1; i++)
		species[i].noneq.assign_Cv_prime(species[i].thermo.Cv[0], species[i].thermo.Cv[1], NONEQ_ROT);

	if (NONEQ_VIB != 0 || NONEQ_E != 0)
	{
		for (i = 0; i <= NS-1; i++)
		{
			if (species[i].basic.type == ELECTRON)
			{
				species[i].noneq.Cv[0]	= 0.0;
			}
		}
	}
}

void FSpecies_updata_Noneq_data2(vector<Species>& species, int NS, int NONEQ_ROT, int NONEQ_VIB, int NONEQ_E)
{
	int	i;

	for (i = 0; i <= NS-1; i++)
		species[i].noneq.assign_Cv_prime(species[i].thermo.Cv[0], species[i].thermo.Cv[1], NONEQ_ROT);

	if (NONEQ_VIB != 0 || NONEQ_E != 0)
	{
		for (i = 0; i <= NS-1; i++)
		{
			if (species[i].basic.type == ELECTRON)
			{
				species[i].noneq.Cv[0]	= 0.0;
			}
		}
	}
}

