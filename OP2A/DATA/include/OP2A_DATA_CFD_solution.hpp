/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 13, 2015
 *      			Author: Minkwan Kim
 *
 * solution_CFD.hpp
 * 			-  
 *  
 */

#ifndef SOLUTION_CFD_HPP_
#define SOLUTION_CFD_HPP_

#include "OP2A_DATA_CFD_setup.hpp"
#include "OP2A_DATA_CFD_mixture.hpp"
#include "OP2A_DATA_CFD_transport_setup.hpp"
#include "OP2A_DATA_CFD_transport_data.hpp"


class SOL_CFD
{
public:
	CFD_variable_setup_ver2				setup;
	CFD_transport_properties_setup_ver2	setup_transport;

	// Cell data
	vector < vector<double> >	Qc;
	vector < vector<double> >	Vc;
	vector < vector<double> >	Wc;
	vector < vector<double> >	Rn;
	vector < vector<double> >	dQ;
	vector < vector<double> >	S_source;
	vector < vector<double> >	Q_new;

	vector < vector<double> >	Qgc;
	vector < vector<double> >	Vgc;
	vector < vector<double> >	Wgc;

	vector < CFD_mixture_data >		mixture_data_c;
	vector < CFD_mixture_data >		mixture_data_gc;

	vector < CFD_transport_properties_data >	transport_prop_c;
	vector < CFD_transport_properties_data >	transport_prop_gc;
	vector < vector < vector<double> > >		Preconditioner;


	// Face data
	vector < vector<double> >	Flux_f_inviscid;
	vector < vector<double> >	Flux_f_viscous;

	// Data for Implicit Method
	vector < vector < vector<double> > >	Jacobian_inviscid_plus;
	vector < vector < vector<double> > >	Jacobian_inviscid_minus;
	vector < vector < vector<double> > >	Jacobian_viscous_plus;
	vector < vector < vector<double> > >	Jacobian_viscous_minus;
	vector < vector < vector<double> > >	Jacobian_source;

	vector < vector < vector<double> > >	dT_dQ;
	vector < vector<double> >				dp_dQ;


	// For Axisymmetric
	vector < double>	div_Vc;


	SOL_CFD();
	~SOL_CFD();


	void allocate_size(int NNM, int NFM, int NCM, int NGM, bool viscous, int time_implicit, bool axisymmetric);
};

#endif /* SOLUTION_CFD_HPP_ */
