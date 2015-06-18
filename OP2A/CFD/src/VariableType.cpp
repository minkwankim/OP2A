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

#include "CFD/include/OP2A_CFD.hpp"

namespace OP2A{
namespace CFD{

void assignVariableType(int NER, int NEV, int NEE)
{
	if (NER == 0)
	{
		if (NEV == 0)
		{
			if (NEE == 0)
			{
				variableType = 1;
			}
			else
			{
				variableType = 2;
			}
		}
		else
		{
			if (NEE == 0)
			{
				variableType = 3;
			}
			else
			{
				variableType = 4;
			}
		}
	}
	else
	{
		if (NEV == 0)
		{
			if (NEE == 0)
			{
				variableType = 5;
			}
			else
			{
				variableType = 6;
			}
		}
		else
		{
			if (NEE == 0)
			{
				variableType = 7;
			}
			else
			{
				variableType = 8;
			}
		}
	}
};

}
}
