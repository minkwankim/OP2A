/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 18, 2015
 *      			Author: Minkwan Kim
 *
 * VariableType.cpp
 * 			-  
 *  
 */

#include "CFD/include/VariableChange.hpp"

namespace OP2A{
namespace CFD{

unsigned int VariableChange::VariableType(int NER, int NEV, int NEE)
{
	unsigned int variabletype = 0;

	if (NER == 0)
	{
		if (NEV == 0)
		{
			if (NEE == 0)
			{
				variabletype = 1;
			}
			else
			{
				variabletype = 2;
			}
		}
		else
		{
			if (NEE == 0)
			{
				variabletype = 3;
			}
			else
			{
				variabletype = 4;
			}
		}
	}
	else
	{
		if (NEV == 0)
		{
			if (NEE == 0)
			{
				variabletype = 5;
			}
			else
			{
				variabletype = 6;
			}
		}
		else
		{
			if (NEE == 0)
			{
				variabletype = 7;
			}
			else
			{
				variabletype = 8;
			}
		}
	}

	return (variabletype);
};



unsigned int VariableChange::NumEnergyModes(unsigned int variableType)
{
	int ne;

	switch (variableType)
	{
	case 1:
		ne = 1;
		break;
	case 2:
		ne = 2;
		break;
	case 3:
		ne = 2;
		break;
	case 4:
		ne = 3;
		break;
	case 5:
		ne = 2;
		break;
	case 6:
		ne = 3;
		break;
	case 7:
		ne = 3;
		break;
	case 8:
		ne = 4;
		break;
	}

	return (ne);
}

}
}
