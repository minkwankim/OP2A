/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 15, 2015
 *      			Author: Minkwan Kim
 *
 * species_info.cpp
 * 			-  
 *  
 */



#include "../include/OP2A_chemistry.hpp"



/*
 * Class NUM 5
 * 	- Species data information
 */
Species_info::Species_info()
{
	NS			= 0;
	electron_ID	= -1;

	// Summary of species info
	numMolecules	= 0;
	for (int s_index = 0; s_index <= MAX_NS_PROBLEM-1; s_index++)
	{
		moleculesID[s_index]	= -1;
	}
}

Species_info::~Species_info()
{

}

void Species_info::assign_data(Species *species, int NS_setup, int showflag)
{
	int		s_index, ns, ns2;
	bool	electron;
	int 	ne, ne2;
	int 	num_molecules, num_molecules2;
	bool 	molecule;


	NS	= NS_setup;

	// Step 1: Obtain molecular information and flagged them as molecules
	// Step 2: Check it has electrons or not
	numMolecules	= 0;
	electron		= false;

	for (s_index = 0; s_index <= NS-1; s_index++)
	{
		if (species[s_index].basic.type == MOLECULE)
		{
			moleculesID[numMolecules]	= s_index;
			numMolecules++;
		}
		else if (species[s_index].basic.type == ELECTRON)
		{
			electron_ID	= s_index;
			electron	= true;
		}
	}

	if (numMolecules >0)	molecule	= true;
	else					molecule = false;

	if (showflag == 1)
	{
		cout << ">>>> Summary of species data <<<<" << endl;
		cout << "================================="	<< endl;
		cout << "Total number of species  : " << NS << endl;
		cout << "Total number of molecules: " << numMolecules << endl;

		if (electron == false)
		{
			cout << "Total number of atoms    : " << NS - numMolecules << endl;
			cout << "  - *it has NO electron"<< endl;
			cout << "=============================" << endl;
		}
		else
		{
			cout << "Total number of atoms    : " << NS - numMolecules -1 << endl;
			cout << "  - *it has electron"<< endl;
			cout << "=============================" << endl;
		}
	}
}


void Species_info::assign_data(vector<Species> species, int NS_setup, int showflag)
{
	int		s_index, ns, ns2;
	bool	electron;
	int 	ne, ne2;
	int 	num_molecules, num_molecules2;
	bool 	molecule;


	NS	= NS_setup;

	// Step 1: Obtain molecular information and flagged them as molecules
	// Step 2: Check it has electrons or not
	numMolecules	= 0;
	electron		= false;

	for (s_index = 0; s_index <= NS-1; s_index++)
	{
		if (species[s_index].basic.type == MOLECULE)
		{
			moleculesID[numMolecules]	= s_index;
			numMolecules++;
		}
		else if (species[s_index].basic.type == ELECTRON)
		{
			electron_ID	= s_index;
			electron	= true;
		}
	}

	if (numMolecules >0)	molecule	= true;
	else					molecule = false;

	if (showflag == 1)
	{
		cout << ">>>> Summary of species data <<<<" << endl;
		cout << "================================="	<< endl;
		cout << "Total number of species  : " << NS << endl;
		cout << "Total number of molecules: " << numMolecules << endl;

		if (electron == false)
		{
			cout << "Total number of atoms    : " << NS - numMolecules << endl;
			cout << "  - *it has NO electron"<< endl;
			cout << "=============================" << endl;
		}
		else
		{
			cout << "Total number of atoms    : " << NS - numMolecules -1 << endl;
			cout << "  - *it has electron"<< endl;
			cout << "=============================" << endl;
		}
	}
}







