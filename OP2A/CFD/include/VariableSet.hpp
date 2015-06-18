/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 16, 2015
 *      			Author: Minkwan Kim
 *
 * VariableSetQ.hpp
 * 			-  
 *  
 */
#ifndef VARIABLESETQ_HPP_
#define VARIABLESETQ_HPP_

#include "CFD/include/OP2A_CFD.hpp"


namespace OP2A{
namespace CFD{

enum FluxType
{
	inv	= 0,
	vis = 1
};

enum FluxCategory
{
	Mass		= 0,
	Momentum	= 1,
	Energy		= 2
};


class CFD_API CFD_VariableSet: public Common::NonInstantiable<CFD_VariableSet>
{
public:
	/*
	 * Basic functions
	 */
	static string var_enerymode(const string& var_name, const CHEM::EnergyMode mode);
	static string var_speciesname(const string& var_name, const string& name);
	static string var_speciesname(const string& var_name, const int& s);
	static string var_stringEnergy(const CHEM::EnergyMode mode);

	static string var_stringQ(const int s, FluxCategory category);
	static string var_stringQ(const string& var_name, FluxCategory category);

	static string var_stringV(const int s, FluxCategory category);
	static string var_stringV(const string& var_name, FluxCategory category);

	static string var_stringW(const int s, FluxCategory category);
	static string var_stringW(const string& var_name, FluxCategory category);

	static string var_stringFluxtype(FluxType type);
	static string var_addUnit(const string& var_name, const string& unit);
	static string var_subTecplot(const string& var_name, const string& sub);
	static string var_subTecplot(const string& var_name, const int& sub);


	/*
	 * Q/V/W variables
	 */
	static string rho_s(const string& name, bool unit);
	static string rho_s(const int& s, bool unit);
	static string velocity(const int& direction, bool unit);
	static string momentum(const int& direction, bool unit);
	static string pressure(bool unit);
	static string T(const CHEM::EnergyMode mode, bool unit);
	static string E(const CHEM::EnergyMode mode, bool unit);



	/*
	 * Flux
	 */
	static string Flux(const int s, FluxCategory category, FluxType type);



	/*
	 * Mixture Data
	 */
	static string rho_mix()	{std::string message = "rho_mix";	return(message);};
	static string R_mix()		{std::string message = "R_mix";		return(message);};
	static string M_mix()		{std::string message = "M_mix";		return(message);};
	static string gamma_mix()	{std::string message = "gamma_mix";	return(message);};
	static string mu_mix()		{std::string message = "mu_mix";	return(message);};
	static string Cv_mix(const CHEM::EnergyMode mode)		{return (var_enerymode("Cv", mode));};
	static string kappa_mix(const CHEM::EnergyMode mode)	{return (var_enerymode("kappa", mode));};


	/*
	 * Mole/Mass fraction
	 */
	static string Xs(const string& name)	{return (var_speciesname("Xs", name));};
	static string Xs(const int& s)			{return (var_speciesname("Xs", s));};

	static string Ys(const string& name)	{return (var_speciesname("Ys", name));};
	static string Ys(const int& s)			{return (var_speciesname("Ys", s));};


	/*
	 * Transport properties
	 */
	static string mu_s(const string& name)	{return (var_speciesname("mu", name));};
	static string mu_s(const int& s)		{return (var_speciesname("mu", s));};

	static string D_s(const string& name)	{return (var_speciesname("D", name));};
	static string D_s(const int& s)		{return (var_speciesname("D", s));};

	static string kappa_s(const CHEM::EnergyMode mode, const string& name)	{return (var_speciesname(var_enerymode("kappa", mode), name));};
	static string kappa_s(const CHEM::EnergyMode mode, const int& s)		{return (var_speciesname(var_enerymode("kappa", mode), s));};


	/*
	 * Derivatives
	 */
	static string dp(const int s, FluxCategory category);
	static string dT(const CHEM::EnergyMode mode, const int s, FluxCategory category);
	static string Jacobian(const int s1, FluxCategory category1, const int s2, FluxCategory category2, FluxType type);

	static Data::DataStorage Q(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage V(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage W(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage MIX(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage Xs(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage Ys(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);

	/*
	 * Utilities
	 */
	static Data::DataStorage FluxTypeVariables(const string& i_VarName, const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage SpeciesTypeVariables(const string& i_VarName, const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage2D JacobianTypeVariables(const string& i_VarName, int mode, const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage2D EnertySpeciesTypeVariables(const string& i_VarName, const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);



	static Data::DataStorage Residue(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage Source(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage dQ(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage Qnew(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage dpdQ(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage DiffusionCoeff(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage ViscosityCoeff(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage DivV(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage FluxInv(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage FluxVis(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);


	static Data::DataStorage2D dTdQ(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage2D dSdQ(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage2D dFinvdQ_plus(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage2D dFvisdQ_plus(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage2D dFinvdQ_minus(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage2D dFvisdQ_minus(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static Data::DataStorage2D thermal_conductivity(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);



	static vector<string> Qstr(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static vector<string> Vstr(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static vector<string> Wstr(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static vector<string> MIXstr(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static vector<string> Xsstr(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
	static vector<string> Ysstr(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis);
};



}
}




#endif /* VARIABLESETQ_HPP_ */
