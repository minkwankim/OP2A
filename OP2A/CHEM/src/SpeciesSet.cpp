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
#include "UTIL/include/OP2A_read_write.hpp"


using namespace std;

namespace OP2A{
namespace CHEM{


SpeciesSet::SpeciesSet(): NS(0), n_atom(0), n_electron(0), n_molecule(0), data_assigned(false), speciesMap(1), NR(0)
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


void SpeciesSet::read_Reaction(const string& file_name)
{
	string text_start 		= "[REACTION START]";
	string text_end 		= "[REACTION END]";
	string text_type 		= "TYPE:";
	string text_kf 			= "FORWARD REACTION COEFFICIENT:";
	string text_Tf 			= "FORWARD REACTION TEMPERATURE:";
	string text_kb 			= "BACKWARD REACTION COEFFICIENT:";
	string text_Tb 			= "BACKWARD REACTION TEMPERATURE:";
	string text_Keq			= "EQUILIBRIUM CONSTANTS:";


	ifstream	reaction_file(file_name.c_str());
	string		line;

	// Step 1: Open file to read
	if(!reaction_file)	throw Common::ExceptionFileSystem (FromHere(), "Could not open species setting file: " + file_name);



	// First reading: finding NR
	NR = 0;
	while (! reaction_file.eof())
	{
		getline(reaction_file, line);
		data_read::read_line(line);

		if (line.compare(0, text_start.size(), text_start) == 0)	// START TO READ REACTION DATA
		{
			while (line.compare(0, text_end.size(), text_end) != 0)
			{
				getline(reaction_file, line);
				data_read::read_line(line);
			}

			NR++;
		}
	}

	reactions.resize(NR);


	// Rewinding for data reading
	reaction_file.clear();
	reaction_file.seekg(0);
	int 	nr = 0;

	string	str_temp;
	int 	temp_int;
	double 	temp_double;
	char   temp_string2[20];
	char   temp_string3[20];
	char   temp_string4[20];


	while (! reaction_file.eof())
	{
		getline(reaction_file, line);
		if (line.compare(0, text_start.size(), text_start) == 0)	// START TO READ REACTION DATA
		{
			while (line.compare(0, text_end.size(), text_end) != 0)
			{
				getline(reaction_file, line);
				data_read::read_line(line);

				temp_int = -1;
				temp_double = -1.0;

				// 1. Reaction name
				str_temp	= data_read::get_data_from_string_string(line, "NAME:");
				if (str_temp !="")
				{
					reactions[nr].name = str_temp;
					reactions[nr].data_processing(str_temp);
				}

				// 2. Type
				data_read::get_data_from_string<int>(line,	"TYPE:", -1, temp_int);
				if (temp_int != -1)
				{
					reactions[nr].type = temp_int;
				}

				// 3. FORWARD REACTION
				if (line.compare(0, text_kf.size(), text_kf) == 0)
				{
					sscanf(line.c_str(),"%s %s %s %lf %lf %lf %lf", temp_string2, temp_string3, temp_string4, &reactions[nr].ForwardReaction.Cf, &reactions[nr].ForwardReaction.nu, &reactions[nr].ForwardReaction.theta, &reactions[nr].ForwardReaction.k);
				}

				// 4. BACKWARD REACTION
				data_read::get_data_from_string<int>(line,	"BACKWARD REACTION RATE METHOD:", -1, temp_int);
				if (temp_int != -1)
				{
					reactions[nr].method_kb = temp_int;
				}

				if (line.compare(0, text_kb.size(), text_kb) == 0)
				{
					sscanf(line.c_str(),"%s	%s %s %lf %lf %lf %lf", temp_string2, temp_string3, temp_string4, &reactions[nr].BackwardReaction.Cf, &reactions[nr].BackwardReaction.nu, &reactions[nr].BackwardReaction.theta, &reactions[nr].BackwardReaction.k);
				}

				// 5. Reference temperature
				data_read::get_data_from_string<double>(line,	"REFERENCE TEMPERATURE:", -1.0, temp_double);
				if (temp_double != -1.0)
				{
					reactions[nr].Tref = temp_double;
				}

				// 6. EQUILIBRIUM CONSTANT DATA
				if (line.compare(0, text_Keq.size(), text_Keq) == 0)
				{
					data_read::get_data_from_string<int>(line,	"EQUILIBRIUM CONSTANTS:", 0, reactions[nr].Keq.model);

					getline(reaction_file, line);	data_read::read_line(line);

					int i	= 0;
					int nn = text_end.size();
					while (line.compare(0, text_end.size(), text_end) != 0)
					{
						if (reactions[nr].Keq.n.size() < (i+1))		reactions[nr].Keq.n.push_back(0.0);
						if (reactions[nr].Keq.A1.size() < (i+1))	reactions[nr].Keq.A1.push_back(0.0);
						if (reactions[nr].Keq.A2.size() < (i+1))	reactions[nr].Keq.A2.push_back(0.0);
						if (reactions[nr].Keq.A3.size() < (i+1))	reactions[nr].Keq.A3.push_back(0.0);
						if (reactions[nr].Keq.A4.size() < (i+1))	reactions[nr].Keq.A4.push_back(0.0);
						if (reactions[nr].Keq.A5.size() < (i+1))	reactions[nr].Keq.A5.push_back(0.0);

						sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf", &reactions[nr].Keq.n[i],
																		&reactions[nr].Keq.A1[i],
																		&reactions[nr].Keq.A2[i],
																		&reactions[nr].Keq.A3[i],
																		&reactions[nr].Keq.A4[i],
																		&reactions[nr].Keq.A5[i]);
						i++;
						getline(reaction_file, line);	data_read::read_line(line);
					}

					if (reactions[nr].Keq.n[i-1] == 0.0)
					{
						reactions[nr].Keq.n.erase(reactions[nr].Keq.n.end());
						reactions[nr].Keq.A1.erase(reactions[nr].Keq.A1.end());
						reactions[nr].Keq.A2.erase(reactions[nr].Keq.A2.end());
						reactions[nr].Keq.A3.erase(reactions[nr].Keq.A3.end());
						reactions[nr].Keq.A4.erase(reactions[nr].Keq.A4.end());
						reactions[nr].Keq.A5.erase(reactions[nr].Keq.A5.end());
					}
				}
			}

			// Update reaction ID
			reactions[nr].ID = nr;
			nr++;
		}
	}
	reaction_file.close();

	// CHECK Read Reaction data
	if (NR != nr)
	{
		std::ostringstream oss;
		oss << "Problem in Reaction DATA: [NR]:" << NR << "  [nr]: " << nr;
		throw Common::ExceptionDimensionMatch (FromHere(), oss.str());
	}

	for (int k = 0; k <= NR-1; k++)
	{
		reactions[k].data_completing(species);
	}
}



void SpeciesSet::showInfoReaction()
{
	cout << ">>>> Summary of Reaction data <<<<" << endl;
	cout << "================================="	<< endl;
	cout << "Total number of Reactions     : " << NR << endl;
	cout << endl;

	int n_diss = 0;
	int n_exch = 0;
	int n_diss_recomb = 0;
	int n_charg_exch = 0;
	int n_e_i_diss = 0;
	int n_e_i_ioni = 0;

	for (int k = 0; k <= NR-1; k++)
	{
		switch(reactions[k].type)
		{
		case ReactionType::DISSOCIATION:
			n_diss++;
			break;

		case ReactionType::EXCHANGE:
			n_exch++;
			break;

		case ReactionType::DISSOCIATIVE_RECOMBINATION:
			n_diss_recomb++;
			break;

		case ReactionType::CHARGE_EXCHANGE:
			n_charg_exch++;
			break;

		case ReactionType::ELECTRON_IMPACT_DISSOCIATION:
			n_e_i_diss++;
			break;

		case ReactionType::ELECTRON_IMPACT_IONIZATION:
			n_e_i_ioni++;
			break;
		}
	}

	cout << "Total number of dissociation \t\t\t: " << n_diss << endl;
	cout << "Total number of exchange \t\t\t: " << n_exch << endl;
	cout << "Total number of dissociative recombination \t: " << n_diss_recomb << endl;
	cout << "Total number of charge exchange \t\t: " << n_charg_exch << endl;
	cout << "Total number of electron impact dissociation \t: " << n_e_i_diss << endl;
	cout << "Total number of electron impact ionization \t: " << n_e_i_ioni << endl;
	cout << "=============================" << endl;
}


}
}
