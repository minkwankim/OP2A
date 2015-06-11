/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 11, 2015
 *      			Author: Minkwan Kim
 *
 * SpeciesTransport.hpp
 * 			-  
 *  
 */
#ifndef SPECIESTRANSPORT_HPP_
#define SPECIESTRANSPORT_HPP_



#include <iostream>
#include <vector>
#include "Common/include/Common.hpp"
#include "CHEM/include/ChemConstants.hpp"



namespace OP2A{
namespace CHEM{


class SpeciesTransport
{
public:

	// Transport properties
	std::vector<double>	Blottner;						// Coefficient for Blottner model: As, Bs, Cs
	std::vector<double>	Sutherland;						// Coefficient for Sutherland model: mu0, T0, S

	double	sigma;								// Collision diameter in Angstrong
	double	T_eps;								// Molecular parameter

	/*
	double	CCS[MAX_NS_PROBLEM][2][4];				// Collision cross section data for neutral collision
	double	OMEGA[MAX_NS_PROBLEM][12][3];			// Collision cross section data for neutral collision
	double	OMEGA_temp[MAX_NS_PROBLEM][12][3];		// Collision cross section data for neutral collision
	*/

	SpeciesTransport();
	~SpeciesTransport();

private:
	bool	data_assignedBlottner;
	bool	data_assignedSutherland;
	bool	data_assignedKinetic;
	bool	data_assignedGupta;
	bool	data_assignedOmega;

	bool	data_assignedTransport;


public:
	void AssignTransport(const std::string& species_name);
	bool is_data_assignedTransport();

private:
	void AssignBlottner(const std::string& species_name);
	void AssignSutherland(const std::string& species_name);
	void AssignKinetic(const std::string& species_name);



};



}
}


#endif /* SPECIESTRANSPORT_HPP_ */
