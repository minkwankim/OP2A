/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 11, 2015
 *      			Author: Minkwan Kim
 *
 * SpeciesTransport.cpp
 * 			-  
 *  
 */




#include <fstream>
#include <sstream>
#include "CHEM/include/ChemConstants.hpp"
#include "CHEM/include/SpeciesTransport.hpp"

#include "Common/include/Exception_FileSystem.hpp"
#include "Common/include/Exception_DimensionMatch.hpp"


using namespace std;

namespace OP2A{
namespace CHEM{



SpeciesTransport::SpeciesTransport():Blottner(3, 0.0), Sutherland(3, 0.0), sigma(0.0), T_eps(0.0),
				data_assignedBlottner(false), data_assignedSutherland(false), data_assignedKinetic(false), data_assignedGupta(false), data_assignedOmega(false), data_assignedTransport(false)
{

}


SpeciesTransport::~SpeciesTransport()
{

}



void SpeciesTransport::AssignBlottner(const std::string& species_name)
{
	string		name_temp;
	ifstream	datafile(OP2A_SPECIES_DATA_BLOT_FILE);

	// Step 1: Open file to read
	if(!datafile)	throw Common::ExceptionFileSystem (FromHere(), "CANNOT FIND A SPECIES DATABASE FILE FOR BLOTTNER MODEL. PLEASE CHECK THE FILE: species_Blottner.dat");


	while (! datafile.eof())
	{
		double		temp1, temp2, temp3;

		datafile >> name_temp	>> temp1 >> temp2 >> temp3;

		int num	= max(species_name.size(), name_temp.size());
		if (name_temp.compare(0, num, species_name)	== 0)
		{
			Blottner[0]				= temp1;
			Blottner[1]				= temp2;
			Blottner[2]				= temp3;
			data_assignedBlottner	= true;
		}
	}

	datafile.close();
}


void SpeciesTransport::AssignSutherland(const std::string& species_name)
{
	string		name_temp;
	ifstream	datafile(OP2A_SPECIES_DATA_SUTH_FILE);

	// Step 1: Open file to read
	if(!datafile)	throw Common::ExceptionFileSystem (FromHere(), "CANNOT FIND A SPECIES DATABASE FILE FOR BLOTTNER MODEL. PLEASE CHECK THE FILE: species_Sutherland.dat");

	while (! datafile.eof())
	{
		double		temp1, temp2, temp3;

		datafile >> name_temp	>> temp1 >> temp2 >> temp3;

		int num	= max(species_name.size(), name_temp.size());
		if (name_temp.compare(0, num, species_name)	== 0)
		{
			Sutherland[0]	= temp1;
			Sutherland[1]	= temp2;
			Sutherland[2]	= temp3;
			data_assignedSutherland	= true;
		}
	}

	datafile.close();
}




void SpeciesTransport::AssignKinetic(const std::string& species_name)
{
	string		name_temp;
	ifstream	datafile(OP2A_SPECIES_DATA_KINE_FILE);

	// Step 1: Open file to read
	if(!datafile)	throw Common::ExceptionFileSystem (FromHere(), "CANNOT FIND A SPECIES DATABASE FILE FOR BLOTTNER MODEL. PLEASE CHECK THE FILE: species_Kinetic.dat");

	while (! datafile.eof())
	{
		double temp1, temp2;

		datafile >> name_temp	>> temp1 >> temp2;

		int num	= max(species_name.size(), name_temp.size());
		if (name_temp.compare(0, num, species_name)	== 0)
		{
			sigma			= temp1;
			T_eps			= temp2;
			data_assignedKinetic	= true;
		}
	}
}


void SpeciesTransport::AssignTransport(const std::string& species_name)
{
	if (data_assignedTransport == false)
	{
		if (data_assignedBlottner == false)		AssignBlottner(species_name);
		if (data_assignedSutherland == false)	AssignSutherland(species_name);
		if (data_assignedKinetic == false)		AssignKinetic(species_name);

		data_assignedTransport = true;
	}
}


bool SpeciesTransport::is_data_assignedTransport()
{
	return (data_assignedTransport);
}


}
}
