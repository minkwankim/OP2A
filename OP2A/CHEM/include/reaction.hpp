/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: May 27, 2014
 *      			Author: Minkwan Kim
 *
 * chemistry.hpp
 * 			-  
 *  
 */

#ifndef _REACTION_HPP_
#define _REACTION_HPP_

#include "OP2A_CHEM_constant.hpp"
#include "species.hpp"
#include "species_data_class.hpp"
#include "species_functions.hpp"

using namespace std;



/*
 * =========================================================
 * 	BASIC Information for reactions
 * 		- Development Status: Version 1.0 (Need to improve)
 * 		- Last modified on: April 7, 2015
 * 						by: Minkwan Kim
 * =========================================================
 */
class REACTION_BASIC
{
public:
	// Reaction name and ID information
	int				ID;					// Reaction ID
	int				type;				// Reaction type
	string			name;


	// REACTANT INFORMATION
	int 			Reactant_num;		// Number of reactant
	vector<int>		Reactant_id;		// ID for global (species_data_entire)
	vector<double>	Reactant_coeff;		// Reactant coefficient, alapha_j


	// PRODUCT INFORMATION
	int				Product_num;		// Number of product
	vector<int>		Product_id;			// ID for global (species_data_entire)
	vector<double>	Product_coeff;		// Product coefficient, beta_j

	REACTION_BASIC();
	~REACTION_BASIC();
};


class REACTION_BASIC_ver2
{
public:
	// Reaction name and ID information
	int				ID;					// Reaction ID
	int				type;				// Reaction type
	string			name;


	// REACTANT INFORMATION
	vector<double>	Reactant_coeff;		// Reactant coefficient, alapha_j

	// PRODUCT INFORMATION
	vector<double>	Product_coeff;		// Product coefficient, beta_j

	REACTION_BASIC_ver2();
	~REACTION_BASIC_ver2();
};




/*
 * =========================================================
 * 	Data for reaction rate
 * 		- Development Status: Version 1.0 (Need to improve)
 * 		- Last modified on: April 7, 2015
 * 						by: Minkwan Kim
 * =========================================================
 */
class REACTION_RATE
{
public:
	int 			method;				// CALCULATION METHOD
	vector<double>	rate_coeff;
	vector<double>	temperature;


	REACTION_RATE();
	~REACTION_RATE();

	double calculate_rate_coefficient(double T);
	double calculate_temperature(vector<double> &T);
};

class REACTION_KEQ
{
public:
	int num_data;
	int model;
	vector<double>	n;
	vector<double>	A1;
	vector<double>	A2;
	vector<double>	A3;
	vector<double>	A4;
	vector<double>	A5;

	REACTION_KEQ();
	~REACTION_KEQ();

	double calculate_logKeq(double T, double n_mix);
	double calculate_dKeq_dT_over_Keq(double T, double n_mix);

};

class REACTION_RATE_COEFF
{
public:
	int 	method;				// CALCULATION METHOD

	double	kf_coeff[4];
	double 	temperature_coeff_f[4];

	double	kb_coeff[4];
	double	temperature_coeff_b[4];

	REACTION_KEQ	Keq_data;
	REACTION_RATE_COEFF();
	~REACTION_RATE_COEFF();
};

void 	calculate_reaction_temperature(REACTION_RATE_COEFF reaction_k, vector<double> &T, double &Tf, double &Tb);
double	cal_kf(REACTION_RATE_COEFF reaction_k, double Tf);
double 	cal_kb(REACTION_RATE_COEFF reaction_k, double n_mix, double Tb);




/*
 * =========================================================
 * 	Data for reaction
 * 		- Development Status: Version 1.0 (Need to improve)
 * 		- Last modified on: April 7, 2015
 * 						by: Minkwan Kim
 * =========================================================
 */

class REACTION : public REACTION_BASIC
{
public:
	double			Tref;
	REACTION_RATE	kf;
	REACTION_RATE	kb;
	REACTION_KEQ	Keq;


	REACTION();
	~REACTION();

	double cal_kf(vector< vector<double> > &T_All, double n_mix);
	double cal_kb(vector< vector<double> > &T_All, double n_mix);
};


class REACTION_SINGLE: public REACTION_BASIC_ver2, public REACTION_RATE_COEFF
{
public:
	double			Tref;
	REACTION_SINGLE();
	~REACTION_SINGLE();
};









/*
 * =================================================
 * 		CHEMICAL REACTION DATA SET
 * =================================================
 */
class	REACTION_DATA
{
public:
	unsigned int	NR;			// Number of total reaction data
	vector<REACTION *>		data_entire;	//

	REACTION_DATA();
	~REACTION_DATA();
};


class	REACTION_DATA_ver2
{
public:
	unsigned int				NR;			// Number of total reaction data
	unsigned int				NS;			// Number of total reaction data

	vector<REACTION_SINGLE>		reaction_k;	//
	vector<SPECIES_ver2> 		species_data;

	vector<vector<double> >		Xf_all;
	vector<vector<double> >		Xf_molecule;
	vector<vector<double> >		Xf_electron;

	vector<vector<double> >		Xb_all;
	vector<vector<double> >		Xb_molecule;
	vector<vector<double> >		Xb_electron;

	REACTION_DATA_ver2();
	~REACTION_DATA_ver2();

	void read_reaction_data(string file_name);
};


void reaction_data_read(string file_name, 	REACTION_DATA &reactions, 		vector<SPECIES_ver2> &species_data);


#endif /* CHEMISTRY_HPP_ */
