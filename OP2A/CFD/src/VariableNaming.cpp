/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 16, 2015
 *      			Author: Minkwan Kim
 *
 * VariableNaming.cpp
 * 			-  
 *  
 */


#include <iostream>     // std::cout, std::ios
#include <sstream>      // std::ostringstream

#include "CFD/include/VariableConstants.hpp"
#include "CFD/include/VariableSet.hpp"
#include "Common/include/Exception_NoSuchValue.hpp"

namespace OP2A{
namespace CFD{

/*
 * Basic functions
 */
string CFD_VariableSet::var_enerymode(const string& var_name, const CHEM::EnergyMode mode)
{
	std::ostringstream message;
	message << var_name << "_" << var_stringEnergy(mode);
	return (message.str());
}


string CFD_VariableSet::var_stringEnergy(const CHEM::EnergyMode mode)
{
	string message;

	switch (mode)
	{
	case CHEM::EnergyMode::E_TRA:
		message = "tra";
		break;

	case CHEM::EnergyMode::E_ROT:
		message = "rot";
		break;

	case CHEM::EnergyMode::E_VIB:
		message = "vib";
		break;
	case CHEM::EnergyMode::E_ELE:
		message  = "e";
		break;

	case CHEM::EnergyMode::E_TR:
		message  = "tr";
		break;

	case CHEM::EnergyMode::E_VE:
		message  = "ve";
		break;

	case CHEM::EnergyMode::E_VEE:
		message = "vee";
		break;
	}

	return (message);
}


string CFD_VariableSet::var_stringQ(const int s, FluxCategory category)
{
	string message;
	switch (category)
	{
	case FluxCategory::Mass:
		message = rho_s(s, false);
		break;

	case FluxCategory::Momentum:
		message = momentum(s, false);
		break;

	case FluxCategory::Energy:
		message = E(static_cast<CHEM::EnergyMode>(s), false);
		break;
	}

	return message;
}


string CFD_VariableSet::var_stringQ(const string& var_name, FluxCategory category)
{
	string message;
	switch (category)
	{
	case FluxCategory::Mass:
		message = rho_s(var_name, false);
		break;
	default:
		throw Common::ExceptionNoSuchValue (FromHere(), "It is not available for Qstring");
		break;
	}

	return message;
}



string CFD_VariableSet::var_stringV(const int s, FluxCategory category)
{
	string message;
	bool use_unit	= true;

	switch (category)
	{
	case FluxCategory::Mass:
		message = rho_s(s, use_unit);
		break;

	case FluxCategory::Momentum:
		message = velocity(s, use_unit);
		break;

	case FluxCategory::Energy:
		message = T(static_cast<CHEM::EnergyMode>(s), use_unit);
		break;
	}

	return message;
}


string CFD_VariableSet::var_stringV(const string& var_name, FluxCategory category)
{
	string message;
	bool use_unit	= true;

	switch (category)
	{
	case FluxCategory::Mass:
		message = rho_s(var_name,  use_unit);
		break;
	default:
		throw Common::ExceptionNoSuchValue (FromHere(), "It is not available for Qstring");
		break;
	}

	return message;
}




string CFD_VariableSet::var_stringW(const int s, FluxCategory category)
{
	string message;
	bool use_unit	= true;

	if (category == FluxCategory::Energy && s == 0)	message	= pressure(use_unit);
	else message = var_stringV(s, category);

	return message;
}


string CFD_VariableSet::var_stringW(const string& var_name, FluxCategory category)
{
	string message = var_stringV(var_name, category);;
	return message;
}






string CFD_VariableSet::var_speciesname(const string& var_name, const string& name)
{
	std::ostringstream message;

	message << var_name << "_" << name;
	return message.str();

}


string CFD_VariableSet::var_speciesname(const string& var_name, const int& s)
{
	std::ostringstream message;

	message << var_name << "_s" << s;
	return message.str();
}


string CFD_VariableSet::var_addUnit(const string& var_name, const string& unit)
{
	std::ostringstream message;

	message << var_name << " [" << unit << "]";
	return message.str();
}


string CFD_VariableSet::var_subTecplot(const string& var_name, const string& sub)
{
	std::ostringstream message;

	message << var_name << "<sub>" << sub << "</sub>";
	return message.str();
}



string CFD_VariableSet::var_subTecplot(const string& var_name, const int& sub)
{
	std::ostringstream message;

	message << var_name << "<sub>" << sub << "</sub>";
	return message.str();
}

string CFD_VariableSet::var_stringFluxtype(FluxType type)
{
	string message;

	switch(type)
	{
	case FluxType::inv:
		message = "inv";
		break;
	case FluxType::vis:
		message = "vis";
		break;
	}

	return message;
}







/*
 * Q/V/W variables
 */

string CFD_VariableSet::rho_s(const string& name, bool unit)
{
	string message;


	if (unit == true)
	{
		message = var_subTecplot("rho", name);
		message = var_addUnit(message, UNIT_DENSITY_SI_TECPLOT);
	}
	else
	{
		message = var_speciesname("rho", name);
	}
	return message;
}

string CFD_VariableSet::rho_s(const int& s, bool unit)
{
	string message;

	message = var_speciesname("rho", s);

	if (unit == true)
		message = var_addUnit(message, UNIT_DENSITY_SI_TECPLOT);

	return message;
}


string CFD_VariableSet::velocity(const int& direction, bool unit)
{
	std::ostringstream message;

	switch (direction)
	{
	case 0:
		message << "u";
		break;

	case 1:
		message << "v";
		break;

	case 2:
		message << "w";
		break;
	}

	if (unit == true) message << "[m/s]";
	return message.str();
}


string CFD_VariableSet::momentum(const int& direction, bool unit)
{
	string message;

	switch (direction)
	{
	case 0:
		message = "rho_u";
		break;

	case 1:
		message = "rho_v";
		break;

	case 2:
		message = "rho_w";
		break;
	}

	if (unit == true)	message = var_addUnit(message, UNIT_MOMENTUM_SI);

	return message;
}



string CFD_VariableSet::pressure(bool unit)
{
	string message = "p";
	if (unit == true)	message = var_addUnit(message, UNIT_PRESSURE_SI);

	return message;
}



string CFD_VariableSet::T(const CHEM::EnergyMode mode, bool unit)
{
	string variablename = "T";
	string modename = var_stringEnergy(mode);
	string message;

	if (unit == true)
	{
		message = var_subTecplot(variablename, modename);
		message = var_addUnit(message, UNIT_TEMPERATURE_SI);
	}
	else
	{
		message = var_enerymode(variablename, mode);
	}

	return message;
}


string CFD_VariableSet::E(const CHEM::EnergyMode mode, bool unit)
{
	string variablename = "E";
	string modename = var_stringEnergy(mode);
	string message;

	if (unit == true)
	{
		message = var_subTecplot(variablename, modename);
		message = var_addUnit(message, UNIT_ENERGY_SI);
	}
	else
	{
		message = var_enerymode(variablename, mode);
	}

	return message;
}




/*
 * Flux
 */
string CFD_VariableSet::Flux(const int s, FluxCategory category, FluxType type)
{
	std::ostringstream message;

	message	<< "F" << var_stringFluxtype(type) << "_";
	message << var_stringQ(s, category);

	return message.str();
}



/*
 * Derivatives
 */
string CFD_VariableSet::dp(const int s, FluxCategory category)
{
	std::ostringstream message;

	message	<< "dp/d" << var_stringQ(s, category);
	return message.str();
}


string CFD_VariableSet::dT(const CHEM::EnergyMode mode, const int s, FluxCategory category)
{
	std::ostringstream message;
	message << "d" << var_enerymode("T", mode) << "/d" << var_stringQ(s, category);
	return message.str();
}


string CFD_VariableSet::Jacobian(const int s1, FluxCategory category1, const int s2, FluxCategory category2, FluxType type)
{
	std::ostringstream message;
	message << "d" << Flux(s1, category1, type) << "/d" << var_stringQ(s2, category2);
	return message.str();
}





}
}
