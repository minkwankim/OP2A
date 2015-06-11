/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 10, 2015
 *      			Author: Minkwan Kim
 *
 * SpeciesBasic.hpp
 * 			-  
 *  
 */
#ifndef SPECIESBASIC_HPP_
#define SPECIESBASIC_HPP_


#include <iostream>
#include <vector>
#include "Common/include/Common.hpp"
#include "CHEM/include/ChemConstants.hpp"



namespace OP2A{
namespace CHEM{



/*
 * Class for Basic Chemistry data
 * @author Minkwan Kim
 * @version 1.0 27/5/2015
 */

enum SpeciesType
{
	Electron	= -1,
	Atom		= 0,
	Molecule	= 1
};




class SpeciesBasic
{
public:

	// Basic information
	std::string 	name;
	int				type;
	double			charge;
	std::string		recombination_to;

	// Basic information
	double	M;									/* ATOMIC MASS, AMU */
	double	m;									/* ATOMIC MASS, kg */
	double	R;									/* Gas constant */
	double	r; 									/* SPECIES RADIUS */

	double	h0;									/* ENTHALPHY OF FORMATION,[J/kg] */
	double	Ds;									/* Dissociation energy [J/kg] */
	double	I;									/* Ionization energy [J/kg] */

	SpeciesBasic();
	explicit SpeciesBasic(const std::string& species_name);

	~SpeciesBasic();

	void AssignDataBasic(const std::string& species_name);
	bool is_assignedBasic();

private:
	bool data_assignedBasic;
};


}
}



#endif /* SPECIESBASIC_HPP_ */
