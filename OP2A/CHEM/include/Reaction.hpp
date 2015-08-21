/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 16, 2015
 *      			Author: Minkwan Kim
 *
 * Resction.hpp
 * 			-  
 *  
 */
#ifndef RESCTION_HPP_
#define RESCTION_HPP_

#include "CHEM/include/Species.hpp"

namespace OP2A{
namespace CHEM{


enum ReactionType
{
	DISSOCIATION					= 0,
	EXCHANGE						= 1,
	DISSOCIATIVE_RECOMBINATION		= 2,
	CHARGE_EXCHANGE					= 3,
	ELECTRON_IMPACT_DISSOCIATION	= 4,
	ELECTRON_IMPACT_IONIZATION		= 5
};

enum KeqModel
{
	PARK85							= 85,
	PARK90							= 90,
	PARK94							= 94,
	LeRC							= 0
};

double Calc_Tc(double T, double Tr, double Tv, double Te, ReactionType type, int direction);


class Rate_coeff_Arrhenius
{
public:
	double Cf;
	double nu;
	double theta;
	double k;

	Rate_coeff_Arrhenius();
	~Rate_coeff_Arrhenius();

	double ReactionRate(double Tc);
};


class EquilibriumConstants
{
public:
	int model;
	std::vector<double>	n;
	std::vector<double>	A1;
	std::vector<double>	A2;
	std::vector<double>	A3;
	std::vector<double>	A4;
	std::vector<double>	A5;

	EquilibriumConstants();
	~EquilibriumConstants();

protected:
	int m_num_data;

public:
	double Keq(double T, double n_mix);

};



class Reaction
{
public:
	// Reaction name and ID information
	int	ID;					// Reaction ID
	int	type;				// Reaction type
	std::string		name;

	// Reactant/product coefficient INFORMATION
	std::vector<double>	alpha;				// Reactant coefficient, alapha_j
	std::vector<double>	beta;				// Product coefficient, beta_j
	std::vector<double>	beta_m_alpha;

	Reaction();
	Reaction(int NS);
	Reaction(std::string data_raw);
	Reaction(std::string data_raw, int i_type);
	Reaction(std::string data_react, std::string data_prod);
	Reaction(std::string data_react, std::string data_prod, int i_type);


	~Reaction();

protected:


public:
	void data_assigning(std::string data_raw);
	void data_processing();
	void data_processing(std::string data_raw);
	void data_completing(const std::vector<Species>& species);


protected:
	std::string m_data_react;
	std::string m_data_prod;

	std::vector<std::string> m_React;
	std::vector<std::string> m_Prod;
	std::vector<double> m_React_c;
	std::vector<double> m_Prod_c;

	bool		m_allocated;
	bool		m_hasdata;

public:
	int method_kb;
	Rate_coeff_Arrhenius	ForwardReaction;
	Rate_coeff_Arrhenius	BackwardReaction;
	EquilibriumConstants	Keq;
	double Tref;

	double kf(double T);
	double kf(double T, double n_mix);

	double kb(double T);
	double kb(double T, double n_mix);
};





}
}


#endif /* RESCTION_HPP_ */
