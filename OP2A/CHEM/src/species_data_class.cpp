/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 13, 2015
 *      			Author: Minkwan Kim
 *
 * species_data_class.cpp
 * 			-  
 *  
 */



#include "../include/OP2A_chemistry.hpp"



/*
 * =================================================
 * 		BASIC SPECIES DATA SET
 * =================================================
 */
SPECIES_DATA_BASIC::SPECIES_DATA_BASIC()
{
	type		= 0;
	NS			= 0;
	data_entire	= NULL;
}

SPECIES_DATA_BASIC::~SPECIES_DATA_BASIC()
{

}


//Internal functions
void SPECIES_DATA_BASIC::allocate_size()
{
	names.resize(NS);		// Name of species
	M.resize(NS);			// Molecular mass in AMU
	m.resize(NS);			// Species mass in kg
	R.resize(NS);			// Species gas constant, [J/kg K]
	h0.resize(NS);			// Species enthalpy of formation
	Cv_t.resize(NS);		// Translational specific heat
	Cv_r.resize(NS);		// Rotational specific heat

	whereis.resize(NS);
}






void SPECIES_DATA_BASIC::allocate_data(vector<SPECIES>	&species_data, int type_flag)
{
	type		= type_flag;
	data_entire	= &species_data;

	switch (type_flag)
	{
	case DATA_ENTIRE:
		NS	= 0;
		for (int s	= 0; s <= data_entire->size()-1; s++)	NS++;
		whereis.resize(NS);

		NS	= 0;
		for (int s	= 0; s <= species_data.size()-1; s++)
		{
			whereis[NS]	= s;
			NS++;
		}
		break;

	case DATA_NEUTRAL:
		NS	= 0;
		for (int s	= 0; s <= species_data.size()-1; s++)
		{
			if (species_data[s].basic_data.q	== 0.0)	NS++;
		}
		whereis.resize(NS);

		NS	= 0;
		for (int s	= 0; s <= species_data.size()-1; s++)
		{
			if (species_data[s].basic_data.q	== 0.0)
			{
				whereis[NS]	= s;
				NS++;
			}
		}
		break;

	case DATA_ION:
		NS	= 0;
		for (int s	= 0; s <= species_data.size()-1; s++)
		{
			if (species_data[s].basic_data.q	!= 0.0	&& species_data[s].basic_data.type != ELECTRON)	NS++;
		}
		whereis.resize(NS);

		NS	= 0;
		for (int s	= 0; s <= species_data.size()-1; s++)
		{
			if (species_data[s].basic_data.q	!= 0.0	&& species_data[s].basic_data.type != ELECTRON)
			{
				whereis[NS]	= s;
				NS++;
			}
		}
		break;

	case DATA_ELECTRON:
		NS	= 0;
		for (int s	= 0; s <= species_data.size()-1; s++)
		{
			if (species_data[s].basic_data.type == ELECTRON)	NS++;
		}
		whereis.resize(NS);

		NS	= 0;
		for (int s	= 0; s <= species_data.size()-1; s++)
		{
			if (species_data[s].basic_data.type == ELECTRON)
			{
				whereis[NS]	= s;
				NS++;
			}
		}
		break;

	case DATA_HEAVYSPECIES:
		NS	= 0;
		for (int s	= 0; s <= species_data.size()-1; s++)
		{
			if (species_data[s].basic_data.type != ELECTRON)	NS++;
		}
		whereis.resize(NS);

		NS	= 0;
		for (int s	= 0; s <= species_data.size()-1; s++)
		{
			if (species_data[s].basic_data.type != ELECTRON)
			{
				whereis[NS]	= s;
				NS++;
			}
		}
		break;


	case DATA_MOLECULES:
		NS	= 0;
		for (int s	= 0; s <= species_data.size()-1; s++)
		{
			if (species_data[s].basic_data.type == MOLECULE)	NS++;
		}
		whereis.resize(NS);

		NS	= 0;
		for (int s	= 0; s <= species_data.size()-1; s++)
		{
			if (species_data[s].basic_data.type == MOLECULE)
			{
				whereis[NS]	= s;
				NS++;
			}
		}
		break;

	case DATA_ATOMS:
		NS	= 0;
		for (int s	= 0; s <= species_data.size()-1; s++)
		{
			if (species_data[s].basic_data.type == ATOM)	NS++;
		}
		whereis.resize(NS);

		NS	= 0;
		for (int s	= 0; s <= species_data.size()-1; s++)
		{
			if (species_data[s].basic_data.type == ATOM)
			{
				whereis[NS]	= s;
				NS++;
			}
		}
		break;
	}


	complete_data();
}






void SPECIES_DATA_BASIC::complete_data()
{
	allocate_size();

	int s_ID;

	for (int s = 0; s <= NS-1; s++)
	{
		s_ID	= whereis[s];

		names[s]	= (*data_entire)[s_ID].basic_data.name;
		M[s]		= (*data_entire)[s_ID].basic_data.M;
		m[s]		= (*data_entire)[s_ID].basic_data.m;
		R[s]		= (*data_entire)[s_ID].thermodynamic_properties.R;
		h0[s]		= (*data_entire)[s_ID].basic_data.h0;
		Cv_t[s]		= (*data_entire)[s_ID].thermodynamic_properties.Cv[0];
		Cv_r[s]		= (*data_entire)[s_ID].thermodynamic_properties.Cv[1];
	}
}

