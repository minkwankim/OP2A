/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 11, 2015
 *      			Author: Minkwan Kim
 *
 * SpeciesSet.cpp
 * 			-  
 *  
 */
#include <fstream>
#include <sstream>
#include <vector>

#include "CHEM/include/SpeciesSet.hpp"
#include "Common/include/Exception_FileSystem.hpp"
#include "Common/include/Exception_DimensionMatch.hpp"


using namespace std;

namespace OP2A{
namespace CHEM{


SpeciesSet::SpeciesSet(): NS(0), n_atom(0), n_electron(0), n_molecule(0), data_assigned(false), speciesMap(1)
{

}

SpeciesSet::SpeciesSet(const std::string& file_name): NS(0), n_atom(0), n_electron(0), n_molecule(0), data_assigned(false), speciesMap(1)
{
	read_SpeciesSet(file_name);
}

SpeciesSet::SpeciesSet(const std::string& file_name, unsigned int ns): NS(ns), n_atom(0), n_electron(0), n_molecule(0), data_assigned(false), speciesMap(ns)
{
	read_SpeciesSet(file_name, ns);
}


SpeciesSet::~SpeciesSet()
{


}




void SpeciesSet::read_SpeciesSet(const std::string& file_name)
{
	ifstream	species_setting(file_name.c_str());

	// Step 1: Open file to read
	if(!species_setting)	throw Common::ExceptionFileSystem (FromHere(), "Could not open species setting file: " + file_name);


	NS	= 0;
	vector<string>	species_name;
	while (! species_setting.eof())
	{
		string name_temp;
		species_setting >> name_temp;

		species_name.push_back(name_temp);
		NS++;
	}
	species_setting.close();

	if (species_name[NS-1] == "")
	{
		NS = NS-1;
		species_name.erase(species_name.end());
	}

	if (species_name.size() != NS)	throw Common::ExceptionDimensionMatch (FromHere(), "Size of species data does not match. Need to check data");

	species.resize(NS);
	speciesMap.reserve(NS);

	for (int s = 0; s <= NS-1; s++)
	{
		species[s].AssignData(species_name[s]);
		speciesMap.insert(species_name[s], s);
	}

	n_atom	= 0;
	n_molecule = 0;
	n_electron = 0;
	for (int s = 0; s <= NS-1; s++)
	{
		switch (species[s].type)
		{
		case SpeciesType::Atom:
			n_atom++;
			break;

		case SpeciesType::Molecule:
			n_molecule++;
			break;

		case SpeciesType::Electron:
			n_electron++;
			break;
		}
	}

	if (n_atom > 0)	atoms.resize(n_atom, NULL);
	if (n_molecule > 0)	molecules.resize(n_molecule, NULL);
	if (n_electron > 0)	electrons.resize(n_electron, NULL);


	n_atom	= 0;
	n_molecule = 0;
	n_electron = 0;
	for (int s = 0; s <= NS-1; s++)
	{
		switch (species[s].type)
		{
		case SpeciesType::Atom:
			atoms[n_atom] = &species[s];
			n_atom++;
			break;

		case SpeciesType::Molecule:
			molecules[n_molecule] = &species[s];
			n_molecule++;
			break;

		case SpeciesType::Electron:
			electrons[n_electron] = &species[s];
			n_electron++;
			break;
		}
	}

	if ((n_atom + n_molecule + n_electron) != NS)	throw Common::ExceptionDimensionMatch (FromHere(), "Size of species data does not match[atom/molecule/electron]. Need to check data");
	else											data_assigned = true;
}



void SpeciesSet::read_SpeciesSet(const std::string& file_name, unsigned int ns)
{
	ifstream	species_setting(file_name.c_str());

	// Step 1: Open file to read
	if(!species_setting)	throw Common::ExceptionFileSystem (FromHere(), "Could not open species setting file: " + file_name);


	NS	= ns;
	vector<string>	species_name(NS);
	species.resize(NS);
	speciesMap.reserve(NS);

	for (int s = 0; s <= NS-1; s++)	species_setting >> species_name[s];
	if (species_name[NS-1] == "")
	{
		NS = NS-1;
		species_name.erase(species_name.end());
	}

	species_setting.close();

	for (int s = 0; s <= NS-1; s++)
	{
		species[s].AssignData(species_name[s]);
		speciesMap.insert(species_name[s], s);
	}



	n_atom	= 0;
	n_molecule = 0;
	n_electron = 0;
	for (int s = 0; s <= NS-1; s++)
	{
		switch (species[s].type)
		{
		case SpeciesType::Atom:
			n_atom++;
			break;

		case SpeciesType::Molecule:
			n_molecule++;
			break;

		case SpeciesType::Electron:
			n_electron++;
			break;
		}
	}


	if (n_atom > 0)	atoms.resize(n_atom, NULL);
	if (n_molecule > 0)	molecules.resize(n_molecule, NULL);
	if (n_electron > 0)	electrons.resize(n_electron, NULL);


	n_atom	= 0;
	n_molecule = 0;
	n_electron = 0;
	for (int s = 0; s <= NS-1; s++)
	{
		switch (species[s].type)
		{
		case SpeciesType::Atom:
			atoms[n_atom] = &species[s];
			n_atom++;
			break;

		case SpeciesType::Molecule:
			molecules[n_molecule] = &species[s];
			n_molecule++;
			break;

		case SpeciesType::Electron:
			electrons[n_electron] = &species[s];
			n_electron++;
			break;
		}
	}

	if ((n_atom + n_molecule + n_electron) != NS)	throw Common::ExceptionDimensionMatch (FromHere(), "Size of species data does not match[atom/molecule/electron]. Need to check data");
	else											data_assigned = true;

}



void SpeciesSet::showInfo()
{
	if (data_assigned != true)	throw Common::ExceptionFileSystem (FromHere(), "Need to read species dataset first!!");

	cout << ">>>> Summary of species data <<<<" << endl;
	cout << "================================="	<< endl;
	cout << "Total number of species     : " << NS << endl;

	cout << "Total number of molecules   : " << n_molecule << endl;
	cout << "	==>List of molecules: ";
	if (n_molecule > 0)
	{
		for (int s = 0; s <= n_molecule-1; s++)
			cout << "[" << molecules[s]->name <<"] ";
		cout << endl;
	}
	else
	{
		cout << "NO Molecules" << endl;
	}

	cout << "Total number of atoms       : " << n_atom << endl;
	cout << "	==>List of atoms    : ";
	if (n_atom > 0)
	{
		for (int s = 0; s <= n_atom-1; s++)
			cout << "[" << atoms[s]->name <<"] ";
		cout << endl;
	}
	else
	{
		cout << "NO Atoms" << endl;
	}

	cout << "Total number of electron    : " << n_electron << endl;
	if (n_electron > 0)	cout << "   **It has electron!" << endl;
	else				cout << "   **It has NO electron!" << endl;

	cout << "=============================" << endl;
}


}
}
