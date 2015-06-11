/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 10, 2015
 *      			Author: Minkwan Kim
 *
 * SpeciesThermo.hpp
 * 			-  
 *  
 */
#ifndef SPECIESTHERMO_HPP_
#define SPECIESTHERMO_HPP_



#include <iostream>
#include <vector>
#include "Common/include/Common.hpp"
#include "CHEM/include/ChemConstants.hpp"



namespace OP2A{
namespace CHEM{

#define MAX_ELETRIC_LEVEL	16
#define	MAX_V_VIB_RELAX		10
#define MAX_VIB_MODE		2



class SpeciesThermo
{
public:
	// THERMODYNAMIC PROPERTIES
	double	Cv_tra;		/* Translational heat capacity */
	double	Cv_rot;		/* Rotational heat capacity */
	double 	Cv_tr;		/* Translational-rotational heat capacity */

	// Rotational Energy
	double	Z_inf;
	double	T_star;
	double	T_ref;
	double	omega;
	double	d_ref;


	// Vibrational Energy
	int 				n_vib_lvl;			/* Number of vibrating mode*/
	std::vector<double>	theta_vib;			/* Characteristic temp of Vibrational energy */

	// Vib-electron relaxation mode
	int					n_vib_ele;
	std::vector<double>	kev1;
	std::vector<double>	kev2;
	std::vector<double>	kev3;
	std::vector<double>	kev4;


	// Electronic Energy
	int 				n_elec_lvl;			/* Number of electronic level */
	std::vector<int> 	g_el;				/* Degeneracies of energy levels */
	std::vector<double>	theta_el;			/* Characteristic electronic temperature */


	// Enthalpy
	int	n_LeRC;
	std::vector< std::vector< double > >	lerc;					/* LeRC curve fits for specific heat, enthalpy and entropy */

	SpeciesThermo();
	~SpeciesThermo();

private:
	bool data_assignedThermo;
	bool data_assignedCv;
	bool data_assignedRot;
	bool data_assignedVib;
	bool data_assignedEle;
	bool data_assignedEnthalpy;
	bool include_kev;


public:	// Internal functions
	void AssignThermoData(const std::string& species_name, const double& R, const int& type);
	bool is_data_assignedThermo();



private:
	void AssignCv(const double& R, const int& type);
	void AssignRot(const std::string& species_name, const int& type);
	void AssignVib(const std::string& species_name, const int& type);
	void AssignEle(const std::string& species_name);
	void AssignEnthalpy(const std::string& species_name);
};



}
}

#endif /* SPECIESTHERMO_HPP_ */
