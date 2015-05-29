/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * CFD_variable_data.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_DATA_CFD_setup.hpp"
#include "../../CFD_new/include/OP2A_CFD_utilities.hpp"

CFD_variable_setup_ver2::CFD_variable_setup_ver2()
{
	NS	= 0;
	ND	= 0;
	NE	= 0;
	VAR	= 0;

	fluid_type	= 0;
	energy_flag	= 0;
	ID_T.resize(4);
	for (int m = 0; m <= 3; m++)	ID_T[m] = 0;
}


CFD_variable_setup_ver2::~CFD_variable_setup_ver2()
{

}



void CFD_variable_setup_ver2::assign_variables(int NER, int NEV, int NEE, int fluid_info)
{
	fluid_type	= fluid_info;

	switch (fluid_type)
	{
	case 0:
		energy_flag = energy_mode_flag(NER, NEV, NEE);
		temperature_mode_table(energy_flag, ID_T);
		NE	= 1 + NER + NEV + NEE;
		break;

	case 1:
		energy_flag = energy_mode_flag(NER, NEV, 0);
		temperature_mode_table(energy_flag, ID_T);
		NE	= 1 + NER + NEV;

		break;

	case -1:
		energy_flag = energy_mode_flag(0, 0, 0);
		temperature_mode_table(energy_flag, ID_T);
		NE	= 1;
		break;
	}

	VAR	= NS + ND + NE;
}
