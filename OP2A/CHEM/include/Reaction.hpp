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
	std::vector<double>	alpha_m_beta;

	Reaction();
	Reaction(int NS);
	Reaction(std::string data_raw);
	Reaction(std::string data_raw, int i_type);
	Reaction(std::string data_react, std::string data_prod);
	Reaction(std::string data_react, std::string data_prod, int i_type);


	~Reaction();

protected:
	void data_assigning(std::string data_raw);
	void data_processing();
	void data_processing(std::string data_raw);

public:
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
};




}
}


#endif /* RESCTION_HPP_ */
