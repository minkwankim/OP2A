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
#include "CHEM/include/Reaction.hpp"
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
	std::vector< Reaction> reactions;

	// Collision data
	std::vector < std::vector<int> > method_collision_integral;	// Collision integral calculation method
															 	// [0]: use CCS data
															 	// [1]: use Wright's DATA(AIAA Journal, Vol 43, No 12, 2005)
															 	// [2]: Attractive force
															 	// [3]: Repulsive force
																// [-1]: NONE

	std::vector< std::vector< std::vector<double> > > 	CCS_A; // Collision cross section data for neutral collision
	std::vector< std::vector< std::vector<double> > > 	CCS_B; // Collision cross section data for neutral collision
	std::vector< std::vector< std::vector<double> > > 	CCS_C; // Collision cross section data for neutral collision
	std::vector< std::vector< std::vector<double> > > 	CCS_D; // Collision cross section data for neutral collision

	std::vector< std::vector< std::vector< std::vector<double> > > > 	Omega; // Collision cross section data for neutral collision
	std::vector< std::vector< std::vector< std::vector<double> > > >	Omega_temp; // Collision cross section data for neutral collision

	SpeciesSet();
	explicit SpeciesSet(const std::string& file_name);
	explicit SpeciesSet(const std::string& file_name, unsigned int ns);


	~SpeciesSet();

private:
	bool data_assigned;

public:
	void read_SpeciesSet(const std::string& file_name);
	void read_SpeciesSet(const std::string& file_name, unsigned int ns);

	void read_Reaction(const std::string& file_name);

	void showInfo();
	void showInfoReaction();


	double pi_Omega0(const int s, const int r, const int l, const double T);
	double pi_Omega1(const int s, const int r, const int l, const double T);
	double pi_Omega2(const int s, const int r, const int l, const double T, const double ne, const double Te);
	double pi_Omega3(const int s, const int r, const int l, const double T, const double ne, const double Te);
	double pi_Omega (const int s, const int r, const int l, const double T, const double ne, const double Te);
	double pi_Omega (const int s, const int r, const int l, const double T);

	double collisionTerm(const int s, const int r, const int l, const double T, const double ne, const double Te);
	double DiffusionBinary_sr(const int s, const int r, const double T, const double ne, const double Te, const double p);

	std::vector<double> Diffusion_s_Gupta(const double T, const double ne, const double Te, const double p, std::vector<double>& Ys);
	double Viscosity_mix_Gupta(const double T, const double ne, const double Te, const double p, std::vector<double>& Ys);
	double ThermalConductivity_tra_Gupta(const double T, const double ne, const double Te, const double p, std::vector<double>& Ys);
	double ThermalConductivity_rot_Gupta(const double T, const double ne, const double Te, const double p, std::vector<double>& Ys);
	double ThermalConductivity_vib_Gupta(const double T, const double ne, const double Te, const double p, std::vector<double>& Ys);
	double ThermalConductivity_e_Gupta(const double T, const double ne, const double Te, const double p, std::vector<double>& Ys);



};




}
}

#endif /* SPECIESSET_HPP_ */
