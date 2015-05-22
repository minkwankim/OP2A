/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Feb 12, 2015
 *      			Author: Minkwan Kim
 *
 * chemisty_model.cpp
 * 			-  
 *  
 */


#include "../include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"


void chemical_production_rate(vector<SPECIES> &species_all, REACTION_DATA &reactions, vector<double> &rho_all, vector< vector<double> > &T_all, vector<double> &n_chem)
{
	int NS	= species_all.size();
	vector <double>	n_all (NS, 0.0);
	vector <double>	rho_M (NS, 0.0);

	for (int s = 0; s <= NS-1; s++)	n_all[s]	= rho_all[s] / species_all[s].basic_data.m;
	for (int s = 0; s <= NS-1; s++)	rho_M[s]	= rho_all[s] / species_all[s].basic_data.M	* 0.001;	// [mol/cm3]


	double n_mix	= 0.0;
	for (int s = 0; s <= NS-1; s++)	n_mix	+= n_all[s];

	for (int s = 0; s <= NS-1; s++)	n_chem[s] = 0.0;

	for (int k = 0; k <= reactions.NR-1; k++)
	{
		double kf	= reactions.data_entire[k]->cal_kf(T_all, n_mix);
		double kb	= reactions.data_entire[k]->cal_kf(T_all, n_mix);

		double Rf	= 1000.0*kf;
		double Rb	= 1000.0*kb;

		for (int s = 0; s <= reactions.data_entire[k]->Reactant_num-1; s++)
		{
			int s_ID	= reactions.data_entire[k]->Reactant_id[s];
			Rf	*= pow(rho_M[s_ID], reactions.data_entire[k]->Reactant_coeff[s]);
		}

		for (int s = 0; s <= reactions.data_entire[k]->Product_num-1; s++)
		{
			int s_ID	= reactions.data_entire[k]->Product_id[s];
			Rb	*= pow(rho_M[s_ID], reactions.data_entire[k]->Product_coeff[s]);
		}


		double Rf_m_Rb	= Rf - Rb;

		// REACTANTS
		for (int s = 0; s <= reactions.data_entire[k]->Reactant_num-1; s++)
		{
			int s_ID		= reactions.data_entire[k]->Reactant_id[s];
			n_chem[s_ID]	-= reactions.data_entire[k]->Reactant_coeff[s]*Rf_m_Rb;
		}

		// PRODUCTS
		for (int s = 0; s <= reactions.data_entire[k]->Product_num-1; s++)
		{
			int s_ID	= reactions.data_entire[k]->Product_id[s];
			n_chem[s_ID]	+= reactions.data_entire[k]->Product_coeff[s]*Rf_m_Rb;
		}
	}


	for (int s = 0; s <= NS-1; s++)	n_chem[s]	/= species_all[s].basic_data.m;

	for (int s = 0; s <= NS-1; s++)
	{
		if (n_chem[s] != n_chem[s]  || fabs(n_chem[s]) == numeric_limits<double>::infinity())
		{
			Error_message_type	error_message;

			error_message.module_name	="CHEMISTRY-REACTION Source Term";
			error_message.location_primary_name	= "N/A";
			error_message.location_secondary_name	= "N/A";
			error_message.message	= " Problem in the calculation of backward reaction rate!";
			error_message.print_message();
		}
	}
}
