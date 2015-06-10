/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 10, 2015
 *      			Author: Minkwan Kim
 *
 * SpeciesThermo.cpp
 * 			-  
 *  
 */



#include <fstream>
#include <sstream>
#include "CHEM/include/SpeciesBasic.hpp"
#include "CHEM/include/SpeciesThermo.hpp"
#include "Common/include/Exception_FileSystem.hpp"
#include "Common/include/Exception_DimensionMatch.hpp"


using namespace std;

namespace OP2A{
namespace CHEM{


SpeciesThermo::SpeciesThermo()
{
	Cv_tra	= 0.0;		/* Translational/Rotational heat capacity */
	Cv_rot	= 0.0;		/* Translational/Rotational heat capacity */
	Cv_tr	= 0.0;

	Z_inf	= 0.0;
	T_star  = 0.0;
	T_ref	= 0.0;
	omega	= 0.0;
	d_ref	= 0.0;


	n_vib_lvl	= 0;
	n_vib_ele	= 0;
	n_elec_lvl	= 0;
	n_LeRC		= 0;

	data_assigned	= false;
	data_assignedCv	= false;
	data_assignedRot	= false;
	data_assignedVib	= false;
	data_assignedEle	= false;
	data_assignedEnthalpy	= false;

	include_kev	= false;
}



SpeciesThermo::~SpeciesThermo()
{

}






/*
 * Internal functions
 */
void SpeciesThermo::AssignCv(const double& R, const int& type)
{
	Cv_tra	= 1.5 * R;
	if (type == SpeciesType::Molecule)	Cv_rot	= R;
	else								Cv_rot	= 0.0;

	Cv_tr	= Cv_rot + Cv_tra;

	data_assignedCv	= true;
}


void SpeciesThermo::AssignRot(const std::string& species_name, const int& type)
{
	int 		num;
	double 		Z_inf_temp, T_star_temp, T_ref_temp, omega_temp, d_ref_temp;


	if (type == SpeciesType::Molecule)
	{
		string		name_temp;
		ifstream	database(OP2A_SPECIES_DATA_ROT_FILE);

		// Step 1: Open file to read
		if(!database)	throw Common::ExceptionFileSystem (FromHere(), "Could not open species database file for Rotational nonequilibrium: species_noneq_rot.dat");

		while (! database.eof())
		{
			database >> name_temp	>> Z_inf_temp >> T_star_temp >> T_ref_temp >> omega_temp >> d_ref_temp;

			num	= max(species_name.size(), name_temp.size());
			if (name_temp.compare(0, num, species_name)	== 0)
			{
				Z_inf	= Z_inf_temp;
				T_star	= T_star_temp;
				T_ref	= T_ref_temp;
				omega	= omega_temp;
				d_ref	= d_ref_temp;

				data_assignedRot	= true;
			}
		}

		database.close();
	}
	else
	{
		data_assignedRot	= true;
	}
}


void SpeciesThermo::AssignVib(const std::string& species_name, const int& type)
{
	if (type == SpeciesType::Molecule)
	{
		string		name_temp;
		double	theta_vib_temp;

		n_vib_lvl	= 0;
		ifstream	database(OP2A_SPECIES_DATA_VIB_FILE);

		// Step 1: Open file to read
		if(!database)	throw Common::ExceptionFileSystem (FromHere(), "Could not open species database file for Vibrational nonequilibrium: species_noneq_vib.dat");

		while (! database.eof())
		{
			database >> name_temp	>> theta_vib_temp;

			int num	= max(species_name.size(), name_temp.size());
			if (name_temp.compare(0, num, species_name)	== 0)
			{
				theta_vib.push_back(theta_vib_temp);
				n_vib_lvl++;
			}
		}

		if (theta_vib.size() == n_vib_lvl)
		{
			data_assignedVib	= true;
		}
		else
		{
			throw Common::ExceptionDimensionMatch (FromHere(), "Size of vibrational non-equilibrium data does not match. Need to check data");
		}
		database.close();



		// Vibrational-electron relaxation
		int		n_kev_temp;
		database.open(OP2A_SPECIES_DATA_K_EV_FILE);
		if(!database)	throw Common::ExceptionFileSystem (FromHere(), "Could not a SPECIES VIBRATION-ELECTRON RELAXATION DATABASE FILE.");

		while (! database.eof())
		{
			database >> name_temp	>> n_kev_temp;
			int num	= max(species_name.size(), name_temp.size());

			if (name_temp.compare(0, num, species_name)	== 0)
			{
				n_vib_ele	= n_kev_temp;
				kev1.resize(n_vib_ele);
				kev2.resize(n_vib_ele);
				kev3.resize(n_vib_ele);
				kev4.resize(n_vib_ele);

				for (int k = 0; k <= n_kev_temp-1; k++)
				{
					database >> kev1[k] >>  kev2[k] >>  kev3[k] >>  kev4[k];
				}

				include_kev = true;
			}
		}
		database.close();

	}
	else
	{
		data_assignedVib	= true;
	}
}


void SpeciesThermo::AssignEle(const std::string& species_name)
{
	ifstream	database(OP2A_SPECIES_DATA_ELE_FILE);

	// Step 1: Open file to read
	if(!database)	throw Common::ExceptionFileSystem (FromHere(), "Could not open an ELECTRONIC STATE DATABASE FILE. PLEASE CHECK THE FILE");


	// Electronic Energy modes
	int		g_i_temp;
	double	theta_el_temp;
	string	name_temp;


	n_elec_lvl	= 0;

	while (! database.eof())
	{
		database >> name_temp	>> theta_el_temp	>> g_i_temp;
		int num	= max(species_name.size(), name_temp.size());

		if (name_temp.compare(0, num, species_name)	== 0)
		{
			theta_el.push_back(theta_el_temp);
			g_el.push_back(g_i_temp);
			n_elec_lvl++;
		}
	}

	if ((theta_el.size() == n_elec_lvl) && (g_el.size() == n_elec_lvl))
	{
		data_assignedEle	= true;
	}
	else
	{
		throw Common::ExceptionDimensionMatch (FromHere(), "Size of electronic state data does not match. Need to check data");
	}
	database.close();
}


void SpeciesThermo::AssignEnthalpy(const std::string& species_name)
{
	ifstream	database(OP2A_SPECIES_DATA_LERC_FILE);

	// Step 1: Open file to read
	if(!database)	throw Common::ExceptionFileSystem (FromHere(), "Could not open a species LeRC DATABASE FILE. PLEASE CHECK THE FILE");


	int 	i_index;
	int		num_lerc	= 0;
	double	data_temp[13];
	string	name_temp;

	n_LeRC	= 0;

	while (! database.eof())
	{
		database >> name_temp	>> data_temp[0] >> data_temp[1] >> data_temp[2] >> data_temp[3] >> data_temp[4] >> data_temp[5] >> data_temp[6] >> data_temp[7] >> data_temp[8] >> data_temp[9] >> data_temp[10] >> data_temp[11] >> data_temp[12];
		int num	= max(species_name.size(), name_temp.size());

		if (name_temp.compare(0, num, species_name)	== 0)
		{
			lerc.push_back(vector<double>(13, 0.0));
			for (int i_index = 0; i_index <= 12; i_index++)
			{
				lerc[n_LeRC][i_index]	= data_temp[i_index];
			}

			n_LeRC++;
		}
	}

	if (lerc.size() == n_LeRC)
	{
		data_assignedEnthalpy	= true;
	}
	else
	{
		throw Common::ExceptionDimensionMatch (FromHere(), "Size of LeRC data does not match. Need to check data");
	}

	database.close();
}




}
}
