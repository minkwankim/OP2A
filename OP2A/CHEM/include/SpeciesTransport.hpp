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
	std::vector<double>	Blottner;							// Coefficient for Blottner model: As, Bs, Cs
	std::vector<double>	Sutherland;							// Coefficient for Sutherland model: mu0, T0, S

	double	sigma;											// Collision diameter in Angstrong
	double	T_eps;											// Molecular parameter


	std::vector<int> method_collision_integral;				 // Collision integral calculation method
															 // [0]: use CCS data
															 // [1]: use Wright's DATA(AIAA Journal, Vol 43, No 12, 2005)
															 // [2]: Attractive force
															 // [3]: Repulsive force

	std::vector< std::vector< std::vector<double> > > 	CCS; // Collision cross section data for neutral collision
	std::vector< std::vector< std::vector<double> > > 	Omega; // Collision cross section data for neutral collision
	std::vector< std::vector< std::vector<double> > > 	Omega_temp; // Collision cross section data for neutral collision

	/*
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

	double piOmega(const int l, const int r, const double T, const double ne, const double Te);
	double piOmega(const int l, const int r, const double T);
	double piOmega0(const int l, const int r, const double T);
	double piOmega1(const int l, const int r, const double T);
	double piOmega2(const int l, const int r, const double T, const double ne, const double Te);
	double piOmega3(const int l, const int r, const double T, const double ne, const double Te);

private:
	void AssignBlottner(const std::string& species_name);
	void AssignSutherland(const std::string& species_name);
	void AssignKinetic(const std::string& species_name);



};



}
}


#endif /* SPECIESTRANSPORT_HPP_ */
