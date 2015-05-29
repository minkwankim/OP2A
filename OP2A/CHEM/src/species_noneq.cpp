/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 15, 2015
 *      			Author: Minkwan Kim
 *
 * species_noneq.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_chemistry.hpp"



/*
 * Class Num 3a
 * 	- Additional variables for Thermal Noneq
 */
Species_NONEQ::Species_NONEQ()
{
	Cv[0]	= 0.0;
	Cv[1]	= 0.0;
}

Species_NONEQ::~Species_NONEQ()
{

}



void Species_NONEQ::assign_Cv_prime(double Cv_tra, double Cv_rot, int NONEQ_ROT)
{
	if (NONEQ_ROT	== 0)
	{
		Cv[0]	= Cv_tra + Cv_rot;
		Cv[1]	= 0.0;
	}
	else
	{
		Cv[0]	= Cv_tra;
		Cv[1]	= Cv_rot;
	}
}



