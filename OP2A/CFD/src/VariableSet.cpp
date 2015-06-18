/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 16, 2015
 *      			Author: Minkwan Kim
 *
 * VariableSet.cpp
 * 			-  
 *  
 */



#include "CFD/include/VariableSet.hpp"
#include "CFD/include/VariableConstants.hpp"

namespace OP2A{
namespace CFD{


/*
 * Basic CFD variables
 */
Data::DataStorage CFD_VariableSet::Q(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int numData;


	numData	= species_set.NS;
	numData	+= ND;
	numData	+= 1;

	if (NER != 0) numData	+= 1;
	if (NEV != 0) numData	+= 1;
	if (NEE != 0) numData	+= 1;

	string variable_name;
	Common::Map1D<std::string, int>	QnameMap (numData);


	/*
	 * Assign Sample data format for
	 * 	 - Q
	 */
		// Species
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		variable_name	= var_stringQ(species_set.species[s].name, FluxCategory::Mass);
		QnameMap.insert(variable_name, s);
	}

	for (int k = 0; k <= ND-1; k++)
	{
		variable_name = var_stringQ(k, FluxCategory::Momentum);
		QnameMap.insert(variable_name, species_set.NS+k);
	}

	int index_ne = species_set.NS+ND;
	variable_name = var_stringQ(CHEM::EnergyMode::E_TRA, FluxCategory::Energy);
	QnameMap.insert(variable_name, index_ne);
	index_ne++;


	// Rotational energy
	if (NER != 0)
	{
		variable_name = var_stringQ(CHEM::EnergyMode::E_ROT, FluxCategory::Energy);
		QnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	// Vibtional energy
	if (NEV != 0)
	{
		variable_name = var_stringQ(CHEM::EnergyMode::E_VE, FluxCategory::Energy);
		QnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	// ELECTRON energy
	if (NEE != 0)
	{
		variable_name = var_stringQ(CHEM::EnergyMode::E_ELE, FluxCategory::Energy);
		QnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	Data::DataStorage	Qsample	(string(NAME_Q), numData, QnameMap);
	return (Qsample);
}

Data::DataStorage CFD_VariableSet::V(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int numData;


	numData	= species_set.NS;
	numData	+= ND;
	numData	+= 1;

	if (NER != 0) numData	+= 1;
	if (NEV != 0) numData	+= 1;
	if (NEE != 0) numData	+= 1;

	string variable_name;
	Common::Map1D<std::string, int>	VnameMap (numData);


	/*
	 * Assign Sample data format for
	 * 	 - Q
	 */
		// Species
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		variable_name	= var_stringV(species_set.species[s].name, FluxCategory::Mass);
		VnameMap.insert(variable_name, s);
	}

	for (int k = 0; k <= ND-1; k++)
	{
		variable_name = var_stringV(k, FluxCategory::Momentum);
		VnameMap.insert(variable_name, species_set.NS+k);
	}

	int index_ne = species_set.NS+ND;
	variable_name = var_stringV(CHEM::EnergyMode::E_TRA, FluxCategory::Energy);
	VnameMap.insert(variable_name, index_ne);
	index_ne++;


	// Rotational energy
	if (NER != 0)
	{
		variable_name = var_stringV(CHEM::EnergyMode::E_ROT, FluxCategory::Energy);
		VnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	// Vibtional energy
	if (NEV != 0)
	{
		variable_name = var_stringV(CHEM::EnergyMode::E_VE, FluxCategory::Energy);
		VnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	// ELECTRON energy
	if (NEE != 0)
	{
		variable_name = var_stringV(CHEM::EnergyMode::E_ELE, FluxCategory::Energy);
		VnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	Data::DataStorage	Vsample	(string(NAME_V), numData, VnameMap);
	return (Vsample);
}

Data::DataStorage CFD_VariableSet::W(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int numData;


	numData	= species_set.NS;
	numData	+= ND;
	numData	+= 1;

	if (NER != 0) numData	+= 1;
	if (NEV != 0) numData	+= 1;
	if (NEE != 0) numData	+= 1;

	string variable_name;
	Common::Map1D<std::string, int>	WnameMap (numData);


	/*
	 * Assign Sample data format for
	 * 	 - Q
	 */
		// Species
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		variable_name	= var_stringW(species_set.species[s].name, FluxCategory::Mass);
		WnameMap.insert(variable_name, s);
	}

	for (int k = 0; k <= ND-1; k++)
	{
		variable_name = var_stringW(k, FluxCategory::Momentum);
		WnameMap.insert(variable_name, species_set.NS+k);
	}

	int index_ne = species_set.NS+ND;
	variable_name = var_stringW(CHEM::EnergyMode::E_TRA, FluxCategory::Energy);
	WnameMap.insert(variable_name, index_ne);
	index_ne++;


	// Rotational energy
	if (NER != 0)
	{
		variable_name = var_stringW(CHEM::EnergyMode::E_ROT, FluxCategory::Energy);
		WnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	// Vibtional energy
	if (NEV != 0)
	{
		variable_name = var_stringW(CHEM::EnergyMode::E_VE, FluxCategory::Energy);
		WnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	// ELECTRON energy
	if (NEE != 0)
	{
		variable_name = var_stringW(CHEM::EnergyMode::E_ELE, FluxCategory::Energy);
		WnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	Data::DataStorage	Wsample	(string(NAME_W), numData, WnameMap);
	return (Wsample);
}

Data::DataStorage CFD_VariableSet::MIX(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int numData;

	int nCv	= 1;
	if (NER != 0)	nCv++;

	int ne	= 1;
	if (NER != 0)	ne++;
	if (NEV != 0)	ne++;
	if (NEE != 0)	ne++;


	numData	= 4 + nCv;
	if(viscous == true)	numData	+= 1 + ne;


	string variable_name;
	Common::Map1D<std::string, int>	MIXnameMap (numData);


	numData = 0;
	variable_name = rho_mix();
	MIXnameMap.insert(variable_name, numData);
	numData++;

	variable_name = R_mix();
	MIXnameMap.insert(variable_name, numData);
	numData++;

	variable_name = M_mix();
	MIXnameMap.insert(variable_name, numData);
	numData++;

	variable_name = gamma_mix();
	MIXnameMap.insert(variable_name, numData);
	numData++;

	if (NER == 0)
	{
		variable_name = Cv_mix(CHEM::EnergyMode::E_TR);
		MIXnameMap.insert(variable_name, numData);
		numData++;
	}
	else
	{
		variable_name = Cv_mix(CHEM::EnergyMode::E_TRA);
		MIXnameMap.insert(variable_name, numData);
		numData++;

		variable_name = Cv_mix(CHEM::EnergyMode::E_ROT);
		MIXnameMap.insert(variable_name, numData);
		numData++;
	}


	if (viscous == true)
	{
		variable_name = mu_mix();
		MIXnameMap.insert(variable_name, numData);
		numData++;

		variable_name = kappa_mix(CHEM::EnergyMode::E_TRA);
		MIXnameMap.insert(variable_name, numData);
		numData++;

		if (NER != 0)
		{
			variable_name = kappa_mix(CHEM::EnergyMode::E_ROT);
			MIXnameMap.insert(variable_name, numData);
			numData++;
		}

		if (NEV != 0)
		{
			variable_name = kappa_mix(CHEM::EnergyMode::E_VIB);
			MIXnameMap.insert(variable_name, numData);
			numData++;
		}

		if (NEE != 0)
		{
			variable_name = kappa_mix(CHEM::EnergyMode::E_ELE);
			MIXnameMap.insert(variable_name, numData);
			numData++;
		}
	}

	Data::DataStorage	MIXsample	(string(NAME_MIX), numData, MIXnameMap);
	return (MIXsample);
}

Data::DataStorage CFD_VariableSet::Xs(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	string variable_name;
	Common::Map1D<std::string, int>	XsnameMap (species_set.NS);

	for (int s = 0; s <= species_set.NS-1; s++)
	{
		variable_name	= Xs(species_set.species[s].name);
		XsnameMap.insert(variable_name, s);
	}

	Data::DataStorage	Xssample	(string(NAME_XS), species_set.NS, XsnameMap);
	return (Xssample);
}

Data::DataStorage CFD_VariableSet::Ys(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	string variable_name;
	Common::Map1D<std::string, int>	YsnameMap (species_set.NS);

	for (int s = 0; s <= species_set.NS-1; s++)
	{
		variable_name	= Ys(species_set.species[s].name);
		YsnameMap.insert(variable_name, s);
	}

	Data::DataStorage	Yssample	(string(NAME_YS), species_set.NS, YsnameMap);
	return (Yssample);
}

Data::DataStorage CFD_VariableSet::Residue(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	string name_front = "Residue_";

	Data::DataStorage o_SampleData = FluxTypeVariables(name_front, species_set, ND, NER, NEV, NEE, viscous, axis);
	o_SampleData.asgName(string(NAME_R));
	return (o_SampleData);
}

Data::DataStorage CFD_VariableSet::Source(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	string name_front = "Source_";

	Data::DataStorage o_SampleData = FluxTypeVariables(name_front, species_set, ND, NER, NEV, NEE, viscous, axis);
	o_SampleData.asgName(string(NAME_S));
	return (o_SampleData);
}

Data::DataStorage CFD_VariableSet::dQ(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	Data::DataStorage	Qsample	= Q(species_set, ND, NER, NEV, NEE, viscous, axis);
	Qsample.asgName(string(NAME_dQ));

	return (Qsample);
}

Data::DataStorage CFD_VariableSet::Qnew(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	Data::DataStorage	Qsample	= Q(species_set, ND, NER, NEV, NEE, viscous, axis);
	Qsample.asgName(string(NAME_Qnew));

	return (Qsample);
}

Data::DataStorage CFD_VariableSet::FluxInv(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	string name_front = "Finv_";

	Data::DataStorage o_SampleData = FluxTypeVariables(name_front, species_set, ND, NER, NEV, NEE, viscous, axis);
	o_SampleData.asgName(string(NAME_FINV));
	return (o_SampleData);
}

Data::DataStorage CFD_VariableSet::FluxVis(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	string name_front = "Fvis_";

	Data::DataStorage o_SampleData = FluxTypeVariables(name_front, species_set, ND, NER, NEV, NEE, viscous, axis);
	o_SampleData.asgName(string(NAME_FVIS));
	return (o_SampleData);
}

Data::DataStorage CFD_VariableSet::dpdQ(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	string name_front = "dp/d";

	Data::DataStorage o_SampleData = FluxTypeVariables(name_front, species_set, ND, NER, NEV, NEE, viscous, axis);
	o_SampleData.asgName(string(NAME_dpdQ));
	return (o_SampleData);
}

Data::DataStorage CFD_VariableSet::DiffusionCoeff(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	string name_front = "D";
	Data::DataStorage	o_SampleData = SpeciesTypeVariables(name_front, species_set, ND, NER, NEV, NEE, viscous, axis);

	o_SampleData.asgName(string(NAME_DS));
	return (o_SampleData);
}

Data::DataStorage CFD_VariableSet::ViscosityCoeff(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	string name_front = "mu";
	Data::DataStorage	o_SampleData = SpeciesTypeVariables(name_front, species_set, ND, NER, NEV, NEE, viscous, axis);

	o_SampleData.asgName(string(NAME_MUS));
	return (o_SampleData);
}

Data::DataStorage CFD_VariableSet::DivV(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	Common::Map1D<std::string, int>	DivVMap (1);

	DivVMap.insert("Divergence V", 0);

	Data::DataStorage sampleData(NAME_DIVV, 1, DivVMap);
	return (sampleData);
}

Data::DataStorage2D CFD_VariableSet::dTdQ(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	Data::DataStorage2D	o_SampleData = JacobianTypeVariables("T", 0, species_set, ND, NER, NEV, NEE, viscous, axis);
	o_SampleData.asgName(NAME_dTdQ);

	return o_SampleData;
}

Data::DataStorage2D CFD_VariableSet::dSdQ(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	Data::DataStorage2D	o_SampleData = JacobianTypeVariables("S", 1, species_set, ND, NER, NEV, NEE, viscous, axis);
	o_SampleData.asgName(NAME_dSdQ);

	return o_SampleData;
}

Data::DataStorage2D CFD_VariableSet::dFinvdQ_plus(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	Data::DataStorage2D	o_SampleData = JacobianTypeVariables("Finv", 1, species_set, ND, NER, NEV, NEE, viscous, axis);
	o_SampleData.asgName(string(NAME_dFINVdQ) + " Plus");

	return o_SampleData;
}

Data::DataStorage2D CFD_VariableSet::dFinvdQ_minus(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	Data::DataStorage2D	o_SampleData = JacobianTypeVariables("Finv", 1, species_set, ND, NER, NEV, NEE, viscous, axis);
	o_SampleData.asgName(string(NAME_dFINVdQ) + " Minus");

	return o_SampleData;
}

Data::DataStorage2D CFD_VariableSet::dFvisdQ_plus(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	Data::DataStorage2D	o_SampleData = JacobianTypeVariables("Fvis", 1, species_set, ND, NER, NEV, NEE, viscous, axis);
	o_SampleData.asgName(string(NAME_dFVISdQ) + " Plus");

	return o_SampleData;
}

Data::DataStorage2D CFD_VariableSet::dFvisdQ_minus(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	Data::DataStorage2D	o_SampleData = JacobianTypeVariables("Fvis", 1, species_set, ND, NER, NEV, NEE, viscous, axis);
	o_SampleData.asgName(string(NAME_dFVISdQ) + " Minus");

	return o_SampleData;
}

Data::DataStorage2D CFD_VariableSet::thermal_conductivity(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	Data::DataStorage2D	o_SampleData = EnertySpeciesTypeVariables("kappa", species_set, ND, NER, NEV, NEE, viscous, axis);

	o_SampleData.asgName(NAME_KAPPAS);
	return o_SampleData;
}




/*
 * Utility functions
 */
Data::DataStorage CFD_VariableSet::FluxTypeVariables(const string& i_VarName, const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int numData;

	numData	= species_set.NS;
	numData	+= ND;
	numData	+= 1;

	if (NER != 0) numData	+= 1;
	if (NEV != 0) numData	+= 1;
	if (NEE != 0) numData	+= 1;


	string variable_name;
	Common::Map1D<std::string, int>	QnameMap (numData);


	/*
	 * Assign Sample data format for
	 * 	 - Q
	 */
		// Species
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		variable_name	= i_VarName + var_stringQ(species_set.species[s].name, FluxCategory::Mass);
		QnameMap.insert(variable_name, s);
	}

	for (int k = 0; k <= ND-1; k++)
	{
		variable_name = i_VarName + var_stringQ(k, FluxCategory::Momentum);
		QnameMap.insert(variable_name, species_set.NS+k);
	}

	int index_ne = species_set.NS+ND;
	variable_name = i_VarName + var_stringQ(CHEM::EnergyMode::E_TRA, FluxCategory::Energy);
	QnameMap.insert(variable_name, index_ne);
	index_ne++;


	// Rotational energy
	if (NER != 0)
	{
		variable_name = i_VarName + var_stringQ(CHEM::EnergyMode::E_ROT, FluxCategory::Energy);
		QnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	// Vibtional energy
	if (NEV != 0)
	{
		variable_name = i_VarName + var_stringQ(CHEM::EnergyMode::E_VE, FluxCategory::Energy);
		QnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	// ELECTRON energy
	if (NEE != 0)
	{
		variable_name = i_VarName + var_stringQ(CHEM::EnergyMode::E_ELE, FluxCategory::Energy);
		QnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	Data::DataStorage	Qsample	("FluxType", numData, QnameMap);
	return (Qsample);
}

Data::DataStorage CFD_VariableSet::SpeciesTypeVariables(const string& i_VarName, const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int numData	= species_set.NS;
	Common::Map1D<std::string, int>	musnameMap (numData);

	for (int s = 0; s <= species_set.NS-1; s++)
	{
		string VarNameTemp;
		VarNameTemp	= var_speciesname(i_VarName, species_set.species[s].name);
		musnameMap.insert(VarNameTemp, s);
	}

	Data::DataStorage	o_SampleData	("SpeciesTypeSample", numData, musnameMap);
	return (o_SampleData);
}

Data::DataStorage2D CFD_VariableSet::JacobianTypeVariables(const string& i_VarName, int mode, const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int ne	= 1;
	if (NER != 0) ne	+= 1;
	if (NEV != 0) ne	+= 1;
	if (NEE != 0) ne	+= 1;


	int numJ = species_set.NS + ND + ne;
	vector<string>	dQname(numJ);

	for (int s = 0; s <= species_set.NS-1; s++)	dQname[s] 					= "d" + var_stringQ(s, FluxCategory::Mass);
	for (int k = 0; k <= ND-1; k++)				dQname[species_set.NS+k]	= "d" + var_stringQ(k, FluxCategory::Momentum);

	ne = 1;
	dQname[species_set.NS+ND] = "d" + var_stringQ(0, FluxCategory::Energy);

	if (NER != 0)
	{
		dQname[species_set.NS+ND+ne] = "d" + var_stringQ(CHEM::EnergyMode::E_ROT, FluxCategory::Energy);
		ne++;
	}

	if (NEV != 0)
	{
		dQname[species_set.NS+ND+ne] = "d" + var_stringQ(CHEM::EnergyMode::E_VIB, FluxCategory::Energy);
		ne++;
	}

	if (NEE != 0)
	{
		dQname[species_set.NS+ND+ne] = "d" + var_stringQ(CHEM::EnergyMode::E_ELE, FluxCategory::Energy);
		ne++;
	}


	int numI;
	switch (mode)
	{
	case 0:
		numI= ne;
		break;
	case 1:
		numI= species_set.NS + ND + ne;
		break;
	}
	vector<string> 	dFname(numI);

	switch (mode)
	{
	case 0:
		ne = 0;
		dFname[ne]	= var_enerymode("d"+i_VarName, CHEM::EnergyMode::E_TRA);	ne++;

		if (NER != 0)
		{
			dFname[ne]	= var_enerymode("d"+i_VarName, CHEM::EnergyMode::E_ROT);
			ne++;
		}

		if (NEV!= 0)
		{
			dFname[ne]	= var_enerymode("d"+i_VarName, CHEM::EnergyMode::E_VIB);
			ne++;
		}

		if (NEE!= 0)
		{
			dFname[ne]	= var_enerymode("d"+i_VarName, CHEM::EnergyMode::E_ELE);
			ne++;
		}
		break;

	case 1:
		for (int s = 0; s <= species_set.NS-1; s++)	dQname[s] 					= "d" + i_VarName + var_stringQ(s, FluxCategory::Mass);
		for (int k = 0; k <= ND-1; k++)				dQname[species_set.NS+k]	= "d" + i_VarName + var_stringQ(k, FluxCategory::Momentum);

		ne = 1;
		dQname[species_set.NS+ND] = "d"+i_VarName + var_stringQ(0, FluxCategory::Energy);

		if (NER != 0)
		{
			dQname[species_set.NS+ND+ne] = "d"+i_VarName + var_stringQ(CHEM::EnergyMode::E_ROT, FluxCategory::Energy);
			ne++;
		}

		if (NEV != 0)
		{
			dQname[species_set.NS+ND+ne] = "d"+i_VarName + var_stringQ(CHEM::EnergyMode::E_VIB, FluxCategory::Energy);
			ne++;
		}

		if (NEE != 0)
		{
			dQname[species_set.NS+ND+ne] = "d"+i_VarName + var_stringQ(CHEM::EnergyMode::E_ELE, FluxCategory::Energy);
			ne++;
		}
		break;
	}


	Common::Map2D<std::string, std::string, int>	o_SampleDataMap (numI*numJ);

	int index = 0;
	for (int j = 0; j <= numJ-1; j++)
	{
		for (int i = 0; i <= numI-1; i++)
		{
			o_SampleDataMap.insert(dFname[i], dQname[j], index);
			index++;
		}
	}


	Data::DataStorage2D	o_SampleData(NAME_dTdQ, numI, numJ, o_SampleDataMap);
	return o_SampleData;
}

Data::DataStorage2D CFD_VariableSet::EnertySpeciesTypeVariables(const string& i_VarName, const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int ne	= 1;
	if (NER != 0) ne	+= 1;
	if (NEV != 0) ne	+= 1;
	if (NEE != 0) ne	+= 1;


	int numJ = species_set.NS;
	vector<string>	dQname(numJ);
	for (int s = 0; s <= species_set.NS-1; s++)	dQname[s] 					= var_stringQ(s, FluxCategory::Mass);


	int numI = ne;
	vector<string> 	dFname(numI);

	ne = 0;
	dFname[ne]	= var_enerymode(i_VarName, CHEM::EnergyMode::E_TRA);	ne++;

	if (NER != 0)
	{
		dFname[ne]	= var_enerymode(i_VarName, CHEM::EnergyMode::E_ROT);
		ne++;
	}

	if (NEV!= 0)
	{
		dFname[ne]	= var_enerymode(i_VarName, CHEM::EnergyMode::E_VIB);
		ne++;
	}

	if (NEE!= 0)
	{
		dFname[ne]	= var_enerymode(i_VarName, CHEM::EnergyMode::E_ELE);
		ne++;
	}


	Common::Map2D<std::string, std::string, int>	o_SampleDataMap (numI*numJ);

	int index = 0;
	for (int j = 0; j <= numJ-1; j++)
	{
		for (int i = 0; i <= numI-1; i++)
		{
			o_SampleDataMap.insert(dFname[i], dQname[j], index);
			index++;
		}
	}

	Data::DataStorage2D	o_SampleData(NAME_dTdQ, numI, numJ, o_SampleDataMap);
	return o_SampleData;
}




/*
 * =============================================
 * Functions for Variable names
 * =============================================
 */

vector<string> CFD_VariableSet::Qstr(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int numData;

	numData	= species_set.NS;
	numData	+= ND;
	numData	+= 1;

	if (NER != 0) numData	+= 1;
	if (NEV != 0) numData	+= 1;
	if (NEE != 0) numData	+= 1;

	vector<string> Qname(numData);


	/*
	 * Assign Sample data format for
	 * 	 - Q
	 */
		// Species
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Qname[s]	= var_stringQ(species_set.species[s].name, FluxCategory::Mass);
	}

	for (int k = 0; k <= ND-1; k++)
	{
		Qname[species_set.NS+k] = var_stringQ(k, FluxCategory::Momentum);
	}

	int index_ne = species_set.NS+ND;
	Qname[index_ne] = var_stringQ(CHEM::EnergyMode::E_TRA, FluxCategory::Energy);
	index_ne++;


	// Rotational energy
	if (NER != 0)
	{
		Qname[index_ne] = var_stringQ(CHEM::EnergyMode::E_ROT, FluxCategory::Energy);
		index_ne++;
	}

	// Vibtional energy
	if (NEV != 0)
	{
		Qname[index_ne] = var_stringQ(CHEM::EnergyMode::E_VE, FluxCategory::Energy);
		index_ne++;
	}

	// ELECTRON energy
	if (NEE != 0)
	{
		Qname[index_ne] = var_stringQ(CHEM::EnergyMode::E_ELE, FluxCategory::Energy);
		index_ne++;
	}

	return (Qname);
}

vector<string> CFD_VariableSet::Vstr(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int numData;

	numData	= species_set.NS;
	numData	+= ND;
	numData	+= 1;

	if (NER != 0) numData	+= 1;
	if (NEV != 0) numData	+= 1;
	if (NEE != 0) numData	+= 1;

	vector<string> Vname(numData);


	/*
	 * Assign Sample data format for
	 * 	 - Q
	 */
		// Species
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Vname[s]	= var_stringV(species_set.species[s].name, FluxCategory::Mass);
	}

	for (int k = 0; k <= ND-1; k++)
	{
		Vname[species_set.NS+k] = var_stringV(k, FluxCategory::Momentum);
	}

	int index_ne = species_set.NS+ND;
	Vname[index_ne] = var_stringV(CHEM::EnergyMode::E_TRA, FluxCategory::Energy);
	index_ne++;


	// Rotational energy
	if (NER != 0)
	{
		Vname[index_ne] = var_stringV(CHEM::EnergyMode::E_ROT, FluxCategory::Energy);
		index_ne++;
	}

	// Vibtional energy
	if (NEV != 0)
	{
		Vname[index_ne] = var_stringV(CHEM::EnergyMode::E_VE, FluxCategory::Energy);
		index_ne++;
	}

	// ELECTRON energy
	if (NEE != 0)
	{
		Vname[index_ne] = var_stringV(CHEM::EnergyMode::E_ELE, FluxCategory::Energy);
		index_ne++;
	}

	return (Vname);
}

vector<string> CFD_VariableSet::Wstr(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int numData;

	numData	= species_set.NS;
	numData	+= ND;
	numData	+= 1;

	if (NER != 0) numData	+= 1;
	if (NEV != 0) numData	+= 1;
	if (NEE != 0) numData	+= 1;

	vector<string> Wname(numData);


	/*
	 * Assign Sample data format for
	 * 	 - Q
	 */
		// Species
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Wname[s]	= var_stringW(species_set.species[s].name, FluxCategory::Mass);
	}

	for (int k = 0; k <= ND-1; k++)
	{
		Wname[species_set.NS+k] = var_stringW(k, FluxCategory::Momentum);
	}

	int index_ne = species_set.NS+ND;
	Wname[index_ne] = var_stringW(CHEM::EnergyMode::E_TRA, FluxCategory::Energy);
	index_ne++;


	// Rotational energy
	if (NER != 0)
	{
		Wname[index_ne] = var_stringW(CHEM::EnergyMode::E_ROT, FluxCategory::Energy);
		index_ne++;
	}

	// Vibtional energy
	if (NEV != 0)
	{
		Wname[index_ne] = var_stringW(CHEM::EnergyMode::E_VE, FluxCategory::Energy);
		index_ne++;
	}

	// ELECTRON energy
	if (NEE != 0)
	{
		Wname[index_ne] = var_stringW(CHEM::EnergyMode::E_ELE, FluxCategory::Energy);
		index_ne++;
	}

	return (Wname);
}

vector<string> CFD_VariableSet::MIXstr(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int numData;

	int nCv	= 1;
	if (NER != 0)	nCv++;

	int ne	= 1;
	if (NER != 0)	ne++;
	if (NEV != 0)	ne++;
	if (NEE != 0)	ne++;


	numData	= 4 + nCv;
	if(viscous == true)	numData	+= 1 + ne;


	vector<string> MIXname(numData);


	numData = 0;
	MIXname[numData] = rho_mix();
	numData++;

	MIXname[numData] = R_mix();
	numData++;

	MIXname[numData] = M_mix();
	numData++;

	MIXname[numData] = gamma_mix();
	numData++;

	if (NER == 0)
	{
		MIXname[numData] = Cv_mix(CHEM::EnergyMode::E_TR);
		numData++;
	}
	else
	{
		MIXname[numData] = Cv_mix(CHEM::EnergyMode::E_TRA);
		numData++;

		MIXname[numData] = Cv_mix(CHEM::EnergyMode::E_ROT);
		numData++;
	}


	if (viscous == true)
	{
		MIXname[numData] = mu_mix();
		numData++;

		MIXname[numData] = kappa_mix(CHEM::EnergyMode::E_TRA);
		numData++;

		if (NER != 0)
		{
			MIXname[numData] = kappa_mix(CHEM::EnergyMode::E_ROT);
			numData++;
		}

		if (NEV != 0)
		{
			MIXname[numData] = kappa_mix(CHEM::EnergyMode::E_VIB);
			numData++;
		}

		if (NEE != 0)
		{
			MIXname[numData] = kappa_mix(CHEM::EnergyMode::E_ELE);
			numData++;
		}
	}

	return (MIXname);
}

vector<string> CFD_VariableSet::Xsstr(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	vector<string> 	Xsname(species_set.NS);

	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Xsname[s]	= Xs(species_set.species[s].name);
	}

	return (Xsname);
}

vector<string> CFD_VariableSet::Ysstr(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	vector<string> 	Ysname(species_set.NS);

	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Ysname[s]	= Ys(species_set.species[s].name);
	}

	return (Ysname);
}


}
}
