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
	int numData;


	numData	= species_set.NS;
	numData	+= ND;
	numData	+= 1;

	if (NER != 0) numData	+= 1;
	if (NEV != 0) numData	+= 1;
	if (NEE != 0) numData	+= 1;

	string name_front = "R_";
	string variable_name;
	Common::Map1D<std::string, int>	QnameMap (numData);


	/*
	 * Assign Sample data format for
	 * 	 - Q
	 */
		// Species
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		variable_name	= name_front + var_stringQ(species_set.species[s].name, FluxCategory::Mass);
		QnameMap.insert(variable_name, s);
	}

	for (int k = 0; k <= ND-1; k++)
	{
		variable_name = name_front + var_stringQ(k, FluxCategory::Momentum);
		QnameMap.insert(variable_name, species_set.NS+k);
	}

	int index_ne = species_set.NS+ND;
	variable_name = name_front + var_stringQ(CHEM::EnergyMode::E_TRA, FluxCategory::Energy);
	QnameMap.insert(variable_name, index_ne);
	index_ne++;


	// Rotational energy
	if (NER != 0)
	{
		variable_name = name_front + var_stringQ(CHEM::EnergyMode::E_ROT, FluxCategory::Energy);
		QnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	// Vibtional energy
	if (NEV != 0)
	{
		variable_name = name_front + var_stringQ(CHEM::EnergyMode::E_VE, FluxCategory::Energy);
		QnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	// ELECTRON energy
	if (NEE != 0)
	{
		variable_name = name_front + var_stringQ(CHEM::EnergyMode::E_ELE, FluxCategory::Energy);
		QnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	Data::DataStorage	Qsample	(string(NAME_R), numData, QnameMap);
	return (Qsample);
}


Data::DataStorage CFD_VariableSet::Source(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int numData;

	numData	= species_set.NS;
	numData	+= ND;
	numData	+= 1;

	if (NER != 0) numData	+= 1;
	if (NEV != 0) numData	+= 1;
	if (NEE != 0) numData	+= 1;

	string name_front = "S_";
	string variable_name;
	Common::Map1D<std::string, int>	QnameMap (numData);


	/*
	 * Assign Sample data format for
	 * 	 - Q
	 */
		// Species
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		variable_name	= name_front + var_stringQ(species_set.species[s].name, FluxCategory::Mass);
		QnameMap.insert(variable_name, s);
	}

	for (int k = 0; k <= ND-1; k++)
	{
		variable_name = name_front + var_stringQ(k, FluxCategory::Momentum);
		QnameMap.insert(variable_name, species_set.NS+k);
	}

	int index_ne = species_set.NS+ND;
	variable_name = name_front + var_stringQ(CHEM::EnergyMode::E_TRA, FluxCategory::Energy);
	QnameMap.insert(variable_name, index_ne);
	index_ne++;


	// Rotational energy
	if (NER != 0)
	{
		variable_name = name_front + var_stringQ(CHEM::EnergyMode::E_ROT, FluxCategory::Energy);
		QnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	// Vibtional energy
	if (NEV != 0)
	{
		variable_name = name_front + var_stringQ(CHEM::EnergyMode::E_VE, FluxCategory::Energy);
		QnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	// ELECTRON energy
	if (NEE != 0)
	{
		variable_name = name_front + var_stringQ(CHEM::EnergyMode::E_ELE, FluxCategory::Energy);
		QnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	Data::DataStorage	Qsample	(string(NAME_S), numData, QnameMap);
	return (Qsample);
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




Data::DataStorage CFD_VariableSet::dpdQ(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int numData;

	numData	= species_set.NS;
	numData	+= ND;
	numData	+= 1;

	if (NER != 0) numData	+= 1;
	if (NEV != 0) numData	+= 1;
	if (NEE != 0) numData	+= 1;

	Common::Map1D<std::string, int>	dpnameMap (numData);


	/*
	 * Assign Sample data format for
	 */
	for (int s = 0; s <= species_set.NS-1; s++)	dpnameMap.insert(dp(s, FluxCategory::Mass), s);
	for (int k = 0; k <= ND-1; k++)				dpnameMap.insert(dp(k, FluxCategory::Momentum), species_set.NS+k);

	string variable_name;
	int index_ne = species_set.NS+ND;

	variable_name = dp(CHEM::EnergyMode::E_TRA, FluxCategory::Energy);
	dpnameMap.insert(variable_name, index_ne);
	index_ne++;


	// Rotational energy
	if (NER != 0)
	{
		variable_name = dp(CHEM::EnergyMode::E_ROT, FluxCategory::Energy);
		dpnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	// Vibtional energy
	if (NEV != 0)
	{
		variable_name = dp(CHEM::EnergyMode::E_VE, FluxCategory::Energy);
		dpnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	// ELECTRON energy
	if (NEE != 0)
	{
		variable_name = dp(CHEM::EnergyMode::E_ELE, FluxCategory::Energy);
		dpnameMap.insert(variable_name, index_ne);
		index_ne++;
	}

	Data::DataStorage	dpsample	(string(NAME_dpdQ), numData, dpnameMap);
	return (dpsample);
}


Data::DataStorage CFD_VariableSet::DiffusionCoeff(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int numData	= species_set.NS;
	Common::Map1D<std::string, int>	DsnameMap (numData);

	for (int s = 0; s <= species_set.NS-1; s++)	DsnameMap.insert(D_s(s), s);

	Data::DataStorage	Dssample	(string(NAME_DS), numData, DsnameMap);
	return (Dssample);
}


Data::DataStorage CFD_VariableSet::ViscosityCoeff(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int numData	= species_set.NS;
	Common::Map1D<std::string, int>	musnameMap (numData);

	for (int s = 0; s <= species_set.NS-1; s++) musnameMap.insert(mu_s(s), s);

	Data::DataStorage	mussample	(string(NAME_MUS), numData, musnameMap);
	return (mussample);
}




Data::DataStorage2D CFD_VariableSet::dTdQ(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int numData = species_set.NS + ND;

	int ne	= 1;
	if (NER != 0) ne	+= 1;
	if (NEV != 0) ne	+= 1;
	if (NEE != 0) ne	+= 1;

	numData += ne;

	Common::Map2D<std::string, std::string, int>	dTdQnameMap (numData*ne);

	vector<string> 	dTname(ne);
	vector<string>	dQname(numData);

	ne = 0;
	dTname[ne]	= var_enerymode("dT", CHEM::EnergyMode::E_TRA);	ne++;
	if (NER != 0)
	{
		dTname[ne]	= var_enerymode("dT", CHEM::EnergyMode::E_ROT);
		ne++;
	}

	if (NEV!= 0)
	{
		dTname[ne]	= var_enerymode("dT", CHEM::EnergyMode::E_VIB);
		ne++;
	}

	if (NEE!= 0)
	{
		dTname[ne]	= var_enerymode("dT", CHEM::EnergyMode::E_ELE);
		ne++;
	}


	numData = 0;
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		dQname[numData] = var_stringQ(s, FluxCategory::Mass);
		numData++;
	}

	for (int k = 0; k <= ND-1; k++)
	{
		dQname[numData] = "d" + var_stringQ(k, FluxCategory::Momentum);
		numData++;
	}

	dQname[numData] = "d" + var_stringQ(0, FluxCategory::Energy);
	numData++;

	if (NER != 0)
	{
		dQname[numData] = "d" + var_stringQ(CHEM::EnergyMode::E_ROT, FluxCategory::Energy);
		numData++;
	}

	if (NEV != 0)
	{
		dQname[numData] = "d" + var_stringQ(CHEM::EnergyMode::E_VIB, FluxCategory::Energy);
		numData++;
	}

	if (NEE != 0)
	{
		dQname[numData] = "d" + var_stringQ(CHEM::EnergyMode::E_ELE, FluxCategory::Energy);
		numData++;
	}


	int index = 0;
	for (int j = 0; j <= numData-1; j++)
	{
		for (int i = 0; i <= ne-1; i++)
		{
			dTdQnameMap.insert(dTname[i], dQname[j], index);
			index++;
		}
	}


	Data::DataStorage2D	dTdQ_sample(NAME_dTdQ, ne, numData, dTdQnameMap);
	return dTdQ_sample;
}



Data::DataStorage2D CFD_VariableSet::thermal_conductivity(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis)
{
	int numData = species_set.NS;

	int ne	= 1;
	if (NER != 0) ne	+= 1;
	if (NEV != 0) ne	+= 1;
	if (NEE != 0) ne	+= 1;


	Common::Map2D<std::string, std::string, int>	kappanameMap (numData*ne);

	vector<string> 	kappa(ne);
	vector<string>	speciesName(numData);

	ne = 0;
	kappa[ne]	= var_enerymode("kappa", CHEM::EnergyMode::E_TRA);	ne++;
	if (NER != 0)
	{
		kappa[ne]	= var_enerymode("kappa", CHEM::EnergyMode::E_ROT);
		ne++;
	}

	if (NEV!= 0)
	{
		kappa[ne]	= var_enerymode("kappa", CHEM::EnergyMode::E_VIB);
		ne++;
	}

	if (NEE!= 0)
	{
		kappa[ne]	= var_enerymode("kappa", CHEM::EnergyMode::E_ELE);
		ne++;
	}


	numData = 0;
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		speciesName[numData] = species_set.species[s].name;
		numData++;
	}

	int index = 0;
	for (int j = 0; j <= numData-1; j++)
	{
		for (int i = 0; i <= ne-1; i++)
		{
			kappanameMap.insert(kappa[ne], speciesName[numData], index);
			index++;
		}
	}


	Data::DataStorage2D	kappa_sample(NAME_KAPPAS, ne, numData, kappanameMap);
	return kappa_sample;
}














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
