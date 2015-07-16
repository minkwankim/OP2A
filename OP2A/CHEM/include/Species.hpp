/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 11, 2015
 *      			Author: Minkwan Kim
 *
 * Species.hpp
 * 			-  
 *  
 */


#ifndef SPECIESOP2A_HPP_
#define SPECIESOP2A_HPP_


#include <iostream>
#include <vector>
#include "Common/include/Common.hpp"

#include "CHEM/include/ChemConstants.hpp"
#include "CHEM/include/SpeciesBasic.hpp"
#include "CHEM/include/SpeciesThermo.hpp"
#include "CHEM/include/SpeciesTransport.hpp"




namespace OP2A{
namespace CHEM{

enum SpeciesViscosityModel
{
	VHS_MODEL			=	0,
	BLOTTNER_MODEL 		=	1,
	SUTHERLAND_MODEL	=	2,
	KINETIC_MODEL		=	3
};

enum SpeciesThermalConductivityModel
{
	EUCKEN_MODEL		=	0
};

enum EnergyMode
{
	E_TRA			=	0,
	E_ROT 			=	1,
	E_VIB			=	2,
	E_ELE			=	3,
	E_TR			=	4,
	E_VE			=	5,
	E_VEE			=	6
};

class Species: public SpeciesBasic, public SpeciesThermo, public SpeciesTransport
{
public:
	Species();
	explicit Species(const std::string& species_name);

	~Species();

	void AssignData(const std::string& species_name);
	bool is_data_assigned();

public:
	double Cv_TRA(const double& T);
	double Cv_ROT(const double& T);
	double Cv_TR(const double& T);

	double Cv_VIB(const double& T);
	double Cv_ELE(const double& T);
	double Cv_VEE(const double& T);
	double Cv_VE(const double& T);

	double e_VIB(const double& T);
	double e_ELE(const double& T);
	double e_VEE(const double& T);
	double e_VE(const double& T);

	double e_internal(const double& Ttra);
	double e_internal(const double& Ttra, const double& Tvib);
	double e_internal(const double& Ttra, const double& Trot, const double& Tvib);
	double e_internal(const double& Ttra, const double& Trot, const double& Tvib, const double& Tele);

	double h(const double& Ttra);
	double h(const double& Ttra, const double& Tvib);
	double h(const double& Ttra, const double& Trot, const double& Tvib);
	double h(const double& Ttra, const double& Trot, const double& Tvib, const double& Tele);

	double viscosityCoeff(const double& T, unsigned int & model);
	double thermalConductivityCoeff(const double& Tvib, const double& Tele, const double& mu,int model,int mode);


private:
	bool data_assigned;

};



}
}

#endif
