/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 16, 2015
 *      			Author: Minkwan Kim
 *
 * CFD_Jacobian_Source_axis.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_CFD_Jacobian.hpp"

void CFD_Jacobian_axis_inviscid(CFD_variable_setup_ver2 setup,	vector<double> &dp_dQ, vector< vector<double> >	&Jacob_axis_inviscid, double S)
{
	int nv	= setup.NS + 1;

	for (int i = 0; i <= setup.VAR-1; i++)	Jacob_axis_inviscid[nv][i]	+= dp_dQ[i]*S;
}


void CFD_Jacobian_axis_viscous(CFD_variable_setup_ver2 setup,	vector<double> &Vc, vector<double> &dp_dQ, CFD_mixture_data &mixture_data, CFD_transport_properties_data &transport_data, vector< vector<double> >	&Jacob_axis_inviscid, double y, double S)
{
	int nv	= setup.NS + 1;

	double lambda 	= -2/3*transport_data.viscosity_coefficient_mixture;
	double temp1	= - (2.0*transport_data.viscosity_coefficient_mixture + lambda) / fabs(y);
	double v_rho	= Vc[nv] / mixture_data.rho;

	for (int s = 0; s <= setup.NS-1; s++)	Jacob_axis_inviscid[nv][s]	+= -temp1*v_rho*S;
	Jacob_axis_inviscid[nv][nv]	+= temp1/mixture_data.rho * S;

	// Need to include further
}
