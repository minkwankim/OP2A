/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 11, 2015
 *      			Author: Minkwan Kim
 *
 * SpeciesSet.hpp
 * 			-  
 *  
 */
#ifndef SPECIESSET_HPP_
#define SPECIESSET_HPP_



#include "CHEM/include/Species.hpp"
#include "Common/include/Map1D.hpp"




namespace OP2A{
namespace CHEM{

class SpeciesSet
{
public:
	unsigned int NS;
	std::vector<Species>	species;
	Common::Map1D<std::string, int>	speciesMap;

	unsigned int n_atom;
	unsigned int n_molecule;
	unsigned int n_electron;

	std::vector<Species*>	atoms;
	std::vector<Species*>	molecules;
	std::vector<Species*>	electrons;

	unsigned int NR;
	std::vector< std::vector <double> > ReactionMatrix;


	SpeciesSet();
	explicit SpeciesSet(const std::string& file_name);
	explicit SpeciesSet(const std::string& file_name, unsigned int ns);


	~SpeciesSet();

private:
	bool data_assigned;

public:
	void read_SpeciesSet(const std::string& file_name);
	void read_SpeciesSet(const std::string& file_name, unsigned int ns);
	void showInfo();

};




}
}

#endif /* SPECIESSET_HPP_ */
