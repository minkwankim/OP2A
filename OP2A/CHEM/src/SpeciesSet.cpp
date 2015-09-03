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
#include "Common/include/MultiDimension.hpp"
#include "Common/include/Exception_FileSystem.hpp"
#include "Common/include/Exception_DimensionMatch.hpp"
#include "UTIL/include/OP2A_read_write.hpp"
#include "Math/include/OP2A_Math.hpp"
#include "Math/include/MathMisc.hpp"

#include "Common/include/Exception_NaNValue.hpp"
#include "Common/include/Exception_NegativeValue.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"

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


	/*
	 * Assign Collision cross section data
	 */
	// 1. Assign Data storage based on NS
	method_collision_integral	= Common::vector_2D<int>(NS, NS, -1);
	CCS_A						= Common::vector_3D<double>(NS, NS, 3, 0.0);
	CCS_B						= Common::vector_3D<double>(NS, NS, 3, 0.0);
	CCS_C						= Common::vector_3D<double>(NS, NS, 3, 0.0);
	CCS_D						= Common::vector_3D<double>(NS, NS, 3, 0.0);
	Omega_temp					= Common::vector_4D<double>(NS, NS, 3, 12, 0.0);
	Omega						= Common::vector_4D<double>(NS, NS, 3, 12, 0.0);


	// 2. Read Collision Cross Section and Omega data
	string		sp1, sp2;
	ifstream	datafile(OP2A_SPECIES_DATA_OMEGA11_FILE);

	// 2.1. Omega 11
	if(!datafile)	throw Common::ExceptionFileSystem (FromHere(), "CANNOT FIND A SPECIES DATABASE FILE FOR COLLISION CROSS SECTION DATA. PLEASE CHECK THE FILE: species_omega11.dat");
	while (! datafile.eof())
	{
		double		data_temp[12];
		double 		data_omega[12];

		int s1 = -1;
		int r1 = -1;

		datafile >> sp1	>> sp2;
		datafile >> data_temp[0]	>> data_temp[1] 	>> data_temp[2] 	>> data_temp[3] 	>> data_temp[4]		>> data_temp[5] 	>> data_temp[6] 	>> data_temp[7] 	>> data_temp[8] 	>> data_temp[9] 	>> data_temp[10] 	>> data_temp[11];
		datafile >> data_omega[0] 	>> data_omega[1] 	>> data_omega[2] 	>> data_omega[3] 	>> data_omega[4]	>> data_omega[5]	>> data_omega[6]	>> data_omega[7]	>> data_omega[8]	>> data_omega[9]	>> data_omega[10]	>> data_omega[11];

		for (int s = 0; s <= NS-1; s++)
		{
			if (sp1 == species[s].name) s1 = s;
			if (sp2 == species[s].name) r1 = s;
		}

		if (s1 != -1 && r1 != -1)
		{
			for (int ii = 0; ii <= 11; ii++)	Omega_temp[s1][r1][1][ii]	= data_temp[ii];
			for (int ii = 0; ii <= 11; ii++)	Omega_temp[r1][s1][1][ii]	= data_temp[ii];

			for (int ii = 0; ii <= 11; ii++)	Omega[s1][r1][1][ii]		= data_omega[ii];
			for (int ii = 0; ii <= 11; ii++)	Omega[r1][s1][1][ii]		= data_omega[ii];

			method_collision_integral[s1][r1]	= 1;
			method_collision_integral[r1][s1]	= 1;
		}
	}
	datafile.close();


	// 2.2. Omega 22
	datafile.open(OP2A_SPECIES_DATA_OMEGA22_FILE);
	if(!datafile)	throw Common::ExceptionFileSystem (FromHere(), "CANNOT FIND A SPECIES DATABASE FILE FOR COLLISION CROSS SECTION DATA. PLEASE CHECK THE FILE: species_omega22.dat");
	while (! datafile.eof())
	{
		double		data_temp[12];
		double 		data_omega[12];

		int s1 = -1;
		int r1 = -1;

		datafile >> sp1	>> sp2;
		datafile >> data_temp[0]	>> data_temp[1] 	>> data_temp[2] 	>> data_temp[3] 	>> data_temp[4]		>> data_temp[5] 	>> data_temp[6] 	>> data_temp[7] 	>> data_temp[8] 	>> data_temp[9] 	>> data_temp[10] 	>> data_temp[11];
		datafile >> data_omega[0] 	>> data_omega[1] 	>> data_omega[2] 	>> data_omega[3] 	>> data_omega[4]	>> data_omega[5]	>> data_omega[6]	>> data_omega[7]	>> data_omega[8]	>> data_omega[9]	>> data_omega[10]	>> data_omega[11];

		for (int s = 0; s <= NS-1; s++)
		{
			if (sp1 == species[s].name) s1 = s;
			if (sp2 == species[s].name) r1 = s;
		}

		if (s1 != -1 && r1 != -1)
		{
			for (int ii = 0; ii <= 11; ii++)	Omega_temp[s1][r1][2][ii]	= data_temp[ii];
			for (int ii = 0; ii <= 11; ii++)	Omega_temp[r1][s1][2][ii]	= data_temp[ii];

			for (int ii = 0; ii <= 11; ii++)	Omega[s1][r1][2][ii]		= data_omega[ii];
			for (int ii = 0; ii <= 11; ii++)	Omega[r1][s1][2][ii]		= data_omega[ii];
		}
	}
	datafile.close();


	// 2.3. CCS: diffusion
	datafile.open(OP2A_SPECIES_DATA_DIFF_FILE);
	if(!datafile)	throw Common::ExceptionFileSystem (FromHere(), "CANNOT FIND A SPECIES DATABASE FILE FOR COLLISION CROSS SECTION DATA. PLEASE CHECK THE FILE: species_diff.dat");
	while (! datafile.eof())
	{
		double		temp1, temp2, temp3, temp4;
		datafile >> sp1	>> sp2 >> temp1 >> temp2 >> temp3 >> temp4;

		int s1 = -1;
		int r1 = -1;

		for (int s = 0; s <= NS-1; s++)
		{
			if (sp1 == species[s].name) s1 = s;
			if (sp2 == species[s].name) r1 = s;
		}

		CCS_A[s1][r1][1]	= temp1;
		CCS_B[s1][r1][1]	= temp2;
		CCS_C[s1][r1][1]	= temp3;
		CCS_D[s1][r1][1]	= temp4;

		CCS_A[r1][s1][1]	= temp1;
		CCS_B[r1][s1][1]	= temp2;
		CCS_C[r1][s1][1]	= temp3;
		CCS_D[r1][s1][1]	= temp4;

		if (method_collision_integral[s1][r1] == -1)	method_collision_integral[s1][r1] = 0;
		if (method_collision_integral[r1][s1] == -1)	method_collision_integral[r1][s1] = 0;
	}
	datafile.close();


	// 2.4. CCS: viscosity
	datafile.open(OP2A_SPECIES_DATA_VIS_FILE);
	if(!datafile)	throw Common::ExceptionFileSystem (FromHere(), "CANNOT FIND A SPECIES DATABASE FILE FOR COLLISION CROSS SECTION DATA. PLEASE CHECK THE FILE: species_vis.dat");
	while (! datafile.eof())
	{
		double		temp1, temp2, temp3, temp4;

		datafile >> sp1	>> sp2 >> temp1 >> temp2 >> temp3 >> temp4;

		int s1 = -1;
		int r1 = -1;

		for (int s = 0; s <= NS-1; s++)
		{
			if (sp1 == species[s].name) s1 = s;
			if (sp2 == species[s].name) r1 = s;
		}

		CCS_A[s1][r1][2]	= temp1;
		CCS_B[s1][r1][2]	= temp2;
		CCS_C[s1][r1][2]	= temp3;
		CCS_D[s1][r1][2]	= temp4;

		CCS_A[r1][s1][2]	= temp1;
		CCS_B[r1][s1][2]	= temp2;
		CCS_C[r1][s1][2]	= temp3;
		CCS_D[r1][s1][2]	= temp4;
	}
	datafile.close();


	// 3. Decide pi_Omega calculation method
	for (int s = 0; s <= NS-1; s++)
	{
		for (int r = 0; r <= NS-1; r++)
		{
			if (method_collision_integral[s][r] == -1)
			{
				if (species[s].charge * species[r].charge < 0.0)
				{
					method_collision_integral[s][r] = 2;
					method_collision_integral[r][s] = 2;
				}

				if (species[s].charge * species[r].charge > 0.0)
				{
					method_collision_integral[s][r] = 3;
					method_collision_integral[r][s] = 3;
				}
			}
		}
	}
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






double SpeciesSet::collisionTerm(const int s, const int r, const int l, const double T, const double ne, const double Te)
{

	double Tp = T;
	if (species[s].type == SpeciesType::Electron || species[r].type == SpeciesType::Electron)	Tp = Te;


	double Rmass;
	Rmass	= 2.0*species[s].M*species[r].M / (species[s].M + species[r].M);

	double aux = sqrt(Rmass/(MATH_PI*Ru*Tp));


	double omega = pi_Omega(s, r, l, T, ne, Te);
	double delta_sr;

	switch (l)
	{
	case 1:
		delta_sr = 8.0/3.0 * aux *omega;
		break;

	case 2:
		delta_sr = 16.0/5.0 * aux *omega;
		break;
	}

	return (delta_sr);
}


double SpeciesSet::DiffusionBinary_sr(const int s, const int r, const double T, const double ne, const double Te, const double p)
{
	double To = T;
	if (species[s].type == SpeciesType::Electron || species[r].type == SpeciesType::Electron)	To = Te;

	double delta;
	double Dsr = 0.0;

	if (s != r)
	{
		delta = collisionTerm(s, r, 1, T, ne, Te);
		Dsr	= C_BOLTZMANN_SI*To / (p*delta);
	}

	return (Dsr);
}



std::vector<double> SpeciesSet::Diffusion_s_Gupta(const double T, const double ne, const double Te, const double p, std::vector<double>& Ys)
{
	vector<double> Ds (NS, 0.0);
	vector<double> gamma_s (NS, 0.0);

	for (int s1 = 0; s1 <= NS-1; s1++)	gamma_s[s1] = Ys[s1]/species[s1].M;

	double gamma_t = 0.0;
	for (int s1 = 0; s1 <= NS-1; s1++)	gamma_t += gamma_s[s1];



	double aux, Dsr;
	for (int s = 0; s <= NS-1; s++)
	{
		aux = 0.0;
		for (int r = 0; r <= NS-1; r++)
		{
			if (r != s)
			{
				Dsr = DiffusionBinary_sr(s, r, T, ne, Te, p);
				aux += gamma_s[r]/Dsr;
			}
		}

		Ds[s] = gamma_t*gamma_t*species[s].M * (1.0 - species[s].M*gamma_s[s]) / aux;


		if (Ds[s] != Ds[s])
		{
			throw Common::ExceptionNaNValue (FromHere(), "NaN value for Diffusion coefficient ");
		}

		if (Ds[s] == numeric_limits<double>::infinity())
		{
			throw Common::ExceptionInfiniteValue (FromHere(), "Infinite Value for Diffusion coefficient");
		}

		if (Ds[s] < 0.0)
		{
			throw Common::ExceptionNegativeValue (FromHere(), "Negative Value for Diffusion coefficient");
		}
	}

	return (Ds);
}




double SpeciesSet::Viscosity_mix_Gupta(const double T, const double ne, const double Te, const double p, std::vector<double>& Ys)
{
	double mu_mix = 0.0;

	vector<double> gamma_s (NS, 0.0);
	for (int s1 = 0; s1 <= NS-1; s1++)	gamma_s[s1] = Ys[s1]/species[s1].M;

	double gamma_t = 0.0;
	for (int s1 = 0; s1 <= NS-1; s1++)	gamma_t += gamma_s[s1];

	double aux;
	for (int s = 0; s <= NS-1; s++)
	{
		aux = 0.0;
		for (int r = 0; r <= NS-1; r++)	aux += gamma_s[r]*	collisionTerm(s, r, 2, T, ne, Te);

		mu_mix += species[s].m*gamma_s[s] / aux;
	}


	if (mu_mix != mu_mix)
	{
		throw Common::ExceptionNaNValue (FromHere(), "NaN value for mixture viscosity coefficient ");
	}

	if (mu_mix == numeric_limits<double>::infinity())
	{
		throw Common::ExceptionInfiniteValue (FromHere(), "Infinite Value for mixture viscosity coefficient");
	}

	if (mu_mix < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Value for mixture viscosity coefficient");
	}

	return (mu_mix);
}


double SpeciesSet::ThermalConductivity_tra_Gupta(const double T, const double ne, const double Te, const double p, std::vector<double>& Ys)
{
	double kappa_tra = 0.0;
	double a_sr;
	double ms_mr;
	double aux1;


	vector<double> gamma_s (NS, 0.0);
	for (int s1 = 0; s1 <= NS-1; s1++)	gamma_s[s1] = Ys[s1]/species[s1].M;

	double gamma_t = 0.0;
	for (int s1 = 0; s1 <= NS-1; s1++)	gamma_t += gamma_s[s1];


	for (int s = 0; s <= NS-1; s++)
	{
		if (species[s].type != SpeciesType::Electron)
		{
			aux1 = 0.0;

			for (int r = 0; r <= NS-1; r++)
			{
				if (species[r].type != SpeciesType::Electron)
				{
					ms_mr = species[s].m / species[r].m;
					a_sr = 1 + (1.0 - ms_mr)*(0.45 - 2.54*ms_mr) / pow(1+ms_mr, 2.0);
				}
				else
				{
					a_sr = 3.54;
				}

				aux1 += gamma_s[r] * collisionTerm(s, r, 2, T, ne, Te);
			}

			kappa_tra += gamma_s[s]/aux1;
		}
	}


	kappa_tra *= 15.0/4.0*C_BOLTZMANN_SI;

	if (kappa_tra != kappa_tra)
	{
		throw Common::ExceptionNaNValue (FromHere(), "NaN value for translational thermal conductivity");
	}

	if (kappa_tra == numeric_limits<double>::infinity())
	{
	throw Common::ExceptionInfiniteValue (FromHere(), "Infinite Value for translational thermal conductivity");
	}

	if (kappa_tra < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Value for translational thermal conductivity");
	}

	return (kappa_tra);
}


double SpeciesSet::ThermalConductivity_rot_Gupta(const double T, const double ne, const double Te, const double p, std::vector<double>& Ys)
{
	double kappa_rot = 0.0;
	double aux1;


	vector<double> gamma_s (NS, 0.0);
	for (int s1 = 0; s1 <= NS-1; s1++)	gamma_s[s1] = Ys[s1]/species[s1].M;

	double gamma_t = 0.0;
	for (int s1 = 0; s1 <= NS-1; s1++)	gamma_t += gamma_s[s1];


	for (int s = 0; s <= NS-1; s++)
	{
		if (species[s].type == SpeciesType::Molecule)
		{
			aux1 = 0.0;

			for (int r = 0; r <= NS-1; r++)
			{
				aux1 += gamma_s[r] * collisionTerm(s, r, 1, T, ne, Te);
			}

			kappa_rot += gamma_s[s]/aux1;
		}
	}

	kappa_rot *= C_BOLTZMANN_SI;

	if (kappa_rot != kappa_rot)
	{
		throw Common::ExceptionNaNValue (FromHere(), "NaN value for rotational thermal conductivity");
	}

	if (kappa_rot == numeric_limits<double>::infinity())
	{
		throw Common::ExceptionInfiniteValue (FromHere(), "Infinite Value for rotational thermal conductivity");
	}

	if (kappa_rot < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Value for rotational thermal conductivity");
	}

	return (kappa_rot);
}



double SpeciesSet::ThermalConductivity_vib_Gupta(const double T, const double ne, const double Te, const double p, std::vector<double>& Ys)
{
	double kappa_vib = 0.0;
	double aux1;
	double aux2;

	vector<double> gamma_s (NS, 0.0);
	for (int s1 = 0; s1 <= NS-1; s1++)	gamma_s[s1] = Ys[s1]/species[s1].M;

	double gamma_t = 0.0;
	for (int s1 = 0; s1 <= NS-1; s1++)	gamma_t += gamma_s[s1];


	for (int s = 0; s <= NS-1; s++)
	{
		if (species[s].type == SpeciesType::Molecule)
		{
			aux1 = 0.0;

			for (int r = 0; r <= NS-1; r++)
			{
				aux1 += gamma_s[r] * collisionTerm(s, r, 1, T, ne, Te);
			}

			aux2 = species[s].Cv_VE(T) / species[s].R;
			kappa_vib += aux2 * (gamma_s[s]/aux1);
		}
	}

	kappa_vib *= C_BOLTZMANN_SI;

	if (kappa_vib != kappa_vib)
	{
		throw Common::ExceptionNaNValue (FromHere(), "NaN value for vibrational thermal conductivity");
	}

	if (kappa_vib == numeric_limits<double>::infinity())
	{
		throw Common::ExceptionInfiniteValue (FromHere(), "Infinite Value for vibrational thermal conductivity");
	}

	if (kappa_vib < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Value for vibrational thermal conductivity");
	}

	return (kappa_vib);
}


double SpeciesSet::ThermalConductivity_e_Gupta(const double T, const double ne, const double Te, const double p, std::vector<double>& Ys)
{
	double kappa_e = 0.0;
	double aux1;

	vector<double> gamma_s (NS, 0.0);
	for (int s1 = 0; s1 <= NS-1; s1++)	gamma_s[s1] = Ys[s1]/species[s1].M;


	int s = -1;
	for (int r = 0; r <= NS-1; r++)
	{
		if (species[r].type == SpeciesType::Electron)
		{
			s = r;
			break;
		}
	}


	aux1 = 0.0;
	for (int r = 0; r <= NS-1; r++)
	{
		aux1 += 1.45 * gamma_s[r] * collisionTerm(s, r, 1, T, ne, Te);
	}

	kappa_e = 15.0/4.0 * C_BOLTZMANN_SI * gamma_s[s] / aux1;

	if (kappa_e != kappa_e)
	{
		throw Common::ExceptionNaNValue (FromHere(), "NaN value for electron thermal conductivity");
	}

	if (kappa_e == numeric_limits<double>::infinity())
	{
		throw Common::ExceptionInfiniteValue (FromHere(), "Infinite Value for electron thermal conductivity");
	}

	if (kappa_e < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Value for electron thermal conductivity");
	}

	return (kappa_e);
}

}
}
