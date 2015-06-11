/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 10, 2015
 *      			Author: Minkwan Kim
 *
 * SpeciesBasic.cpp
 * 			-  
 *  
 */

#include <fstream>
#include <sstream>
#include "CHEM/include/SpeciesBasic.hpp"
#include "Common/include/Exception_FileSystem.hpp"

using namespace std;

namespace OP2A{
namespace CHEM{




SpeciesBasic::SpeciesBasic()
{
	type 	= -2;
	charge	= 0;

	M		= 0.0;									/* ATOMIC MASS, AMU */
	m		= 0.0;									/* ATOMIC MASS, kg */
	R		= 0.0;									/* Gas constant */
	r		= 0.0; 									/* SPECIES RADIUS */

	h0		= 0.0;									/* ENTHALPHY OF FORMATION,[J/kg] */
	Ds		= 0.0;									/* Dissociation energy [J/kg] */
	I		= 0.0;									/* Ionization energy [J/kg] */

	data_assignedBasic = false;
}


SpeciesBasic::SpeciesBasic(const string& species_name)
{
	AssignDataBasic(species_name);
}


SpeciesBasic::~SpeciesBasic()
{

}


/*
 * Internal Functions
 */
void SpeciesBasic::AssignDataBasic(const string& species_name)
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


	string		line;
	ifstream	database(OP2A_SPECIES_DATA_BASE_FILE);

	name	= species_name;

	// Step 1: Open file to read
	if(!database)	throw Common::ExceptionFileSystem (FromHere(), "Could not open species database file");

	while (! database.eof())
	{
		database >> name_temp	>> M_temp	>> h0_temp	>> Ds_temp	>> I_temp	>> q_temp	>> type_temp	>> recombination_temp;

		int num	= max(name.size(), name_temp.size());
		if (name_temp.compare(0, num, name)	== 0)
		{
			M					= M_temp;
			m					= M * AMU_SI;
			h0					= h0_temp;
			Ds					= Ds_temp;
			I					= I_temp;
			charge				= q_temp;
			type				= type_temp;
			recombination_to	= recombination_temp;
			R					= Ru / M;

			data_assignedBasic		= true;
		}
	}

	database.close();
}


bool SpeciesBasic::is_assignedBasic()
{
	return (data_assignedBasic);
}



}
}
