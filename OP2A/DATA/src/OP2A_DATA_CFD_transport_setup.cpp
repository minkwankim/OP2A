/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_DATA_CFD_transport_setup.cpp
 * 			-  
 *  
 */

#ifndef OP2A_DATA_CFD_TRANSPORT_SETUP_CPP_
#define OP2A_DATA_CFD_TRANSPORT_SETUP_CPP_



#include "../include/OP2A_DATA_CFD_transport_setup.hpp"



CFD_transport_properties_setup_ver2::CFD_transport_properties_setup_ver2()
{
	viscosity_model		= -1;
	conductivity_model	= -1;
	mixing_rule			= -1;
	Le					= 0.0;
}

CFD_transport_properties_setup_ver2::~CFD_transport_properties_setup_ver2()
{

}


#endif /* OP2A_DATA_CFD_TRANSPORT_SETUP_CPP_ */
