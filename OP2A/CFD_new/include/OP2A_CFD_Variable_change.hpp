/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_CFD_Variable_change.hpp
 * 			-  
 *  
 */

#ifndef OP2A_CFD_VARIABLE_CHANGE_HPP_
#define OP2A_CFD_VARIABLE_CHANGE_HPP_



#include "../../DATA/include/OP2A_DATA_CFD_solution.hpp"
#include "../../CHEM/include/OP2A_chemistry.hpp"


class fn_Evib
{
public:
	SPECIES_DATA_BASIC			species;
	vector<double>				rhos;
	int 						NS;
	double						Eve;

	fn_Evib();
	~fn_Evib();

	void	create(SPECIES_DATA_BASIC	&species_data, vector<double> &rho, int ns, double E);
	double	f(double T);
	double	df(double T);
};

class fn_Evib2
{
public:
	SPECIES_DATA_BASIC			species;
	vector<double>				rhos;
	int 						NS;
	double						Eve;

	fn_Evib2();
	~fn_Evib2();

	void	create(SPECIES_DATA_BASIC	&species_data, vector<double> &rho, int ns, double E);
	double	f(double T);
	double	df(double T);
};



void CFD_V_to_Q(vector<double> V, vector<double> &Q, CFD_variable_setup_ver2 setup, CFD_mixture_data variable_data, SPECIES_DATA_BASIC &species);
void CFD_Q_to_V(vector<double> Q, vector<double> &V, CFD_variable_setup_ver2 setup, CFD_mixture_data variable_data, SPECIES_DATA_BASIC &species);
void CFD_V_to_W(vector<double> V, vector<double> &W, CFD_variable_setup_ver2 setup, CFD_mixture_data variable_data, SPECIES_DATA_BASIC &species);
void CFD_W_to_V(vector<double> W, vector<double> &V, CFD_variable_setup_ver2 setup, CFD_mixture_data variable_data, SPECIES_DATA_BASIC &species);
void CFD_Calculate_all_required_variables_inviscid(SOL_CFD &Solution, int NCM, SPECIES_DATA_BASIC	&species, int NT, int mode);



#endif /* OP2A_CFD_VARIABLE_CHANGE_HPP_ */
