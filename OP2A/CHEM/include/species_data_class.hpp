/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 13, 2015
 *      			Author: Minkwan Kim
 *
 * species_data_class.hpp
 * 			-  
 *  
 */

#ifndef SPECIES_DATA_CLASS_HPP_
#define SPECIES_DATA_CLASS_HPP_

#include <vector>
#include <limits>


#include "species.hpp"
#include "species_functions.hpp"



/*
 * =================================================
 * 		BASIC SPECIES DATA SET
 * =================================================
 */
class	SPECIES_DATA_BASIC
{
public:
	int	type;
	unsigned int	NS;			// Number of species

	vector<string>	names;		// Name of species
	vector<double>	M;			// Molecular mass in AMU
	vector<double>	m;			// Species mass in kg
	vector<double>	R;			// Species gas constant, [J/kg K]
	vector<double>	h0;			// Species enthalpy of formation
	vector<double>	Cv_t;		// Translational specific heat
	vector<double>	Cv_r;		// Rotational specific heat



	vector<int>			whereis;
	vector<SPECIES>	*	data_entire;	//

	SPECIES_DATA_BASIC();
	~SPECIES_DATA_BASIC();

	//Internal functions
	void allocate_size();
	void allocate_data(vector<SPECIES>	&species_data_entire, int type_flag);
	void complete_data();

};









#endif /* SPECIES_DATA_CLASS_HPP_ */
