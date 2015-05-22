/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_utilities_energymode.cpp
 * 			-  
 *  
 */


#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include "../../CHEM/include/OP2A_CHEM_constant.hpp"
using namespace std;


int	energy_mode_flag(int NER, int NEV, int NEE)
{
	int flag;

	if (NER == 0)
	{
		if (NEV == 0)
		{
			if (NEE == 0)	flag = 0;
			else			flag = 1;
		}
		else
		{
			if (NEE == 0)	flag = 2;
			else			flag = 3;
		}
	}
	else
	{
		if (NEV == 0)
		{
			if (NEE == 0)	flag = 4;
			else			flag = 5;
		}
		else
		{
			if (NEE == 0)	flag = 6;
			else			flag = 7;
		}
	}

	return (flag);
}


void temperature_mode_table(int flag, vector<int> &ID_T)
{
	ID_T.resize(4);

	ID_T[TRA]	= 0;
	switch(flag)
	{
	case 0:
		ID_T[ROT]	= 0;
		ID_T[VIB]	= 0;
		ID_T[ELE]	= 0;
		break;

	case 1:
		ID_T[ROT]	= 0;
		ID_T[VIB]	= 0;
		ID_T[ELE]	= 1;
		break;

	case 2:
		ID_T[ROT]	= 0;
		ID_T[VIB]	= 1;
		ID_T[ELE]	= 0;
		break;

	case 3:
		ID_T[ROT]	= 0;
		ID_T[VIB]	= 1;
		ID_T[ELE]	= 2;
		break;

	case 4:
		ID_T[ROT]	= 1;
		ID_T[VIB]	= 0;
		ID_T[ELE]	= 0;
		break;

	case 5:
		ID_T[ROT]	= 1;
		ID_T[VIB]	= 0;
		ID_T[ELE]	= 2;
		break;

	case 6:
		ID_T[ROT]	= 1;
		ID_T[VIB]	= 2;
		ID_T[ELE]	= 0;
		break;

	case 7:
		ID_T[ROT]	= 1;
		ID_T[VIB]	= 1;
		ID_T[ELE]	= 1;
		break;
	}
}
