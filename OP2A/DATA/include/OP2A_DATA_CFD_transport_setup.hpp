/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_DATA_CFD_transport_setup.hpp
 * 			-  
 *  
 */

#ifndef OP2A_DATA_CFD_TRANSPORT_SETUP_HPP_
#define OP2A_DATA_CFD_TRANSPORT_SETUP_HPP_


class CFD_transport_properties_setup_ver2
{
public:
	int 	viscosity_model;
	int 	conductivity_model;
	int 	mixing_rule;
	double	Le;

	CFD_transport_properties_setup_ver2();
	~CFD_transport_properties_setup_ver2();
};


#endif /* OP2A_DATA_CFD_TRANSPORT_SETUP_HPP_ */
