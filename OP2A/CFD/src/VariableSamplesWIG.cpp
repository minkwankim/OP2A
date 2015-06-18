/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 17, 2015
 *      			Author: Minkwan Kim
 *
 * VariableSamples.cpp
 * 			-  
 *  
 */


#include "CFD/include/VariableSamplesWIG.hpp"
#include "CFD/include/VariableConstants.hpp"


namespace OP2A{
namespace CFD{



Data::DataStorageVector<Data::DataStorage> CFD_DataTemplateWIG::CellData1D(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis, int time_integration)
{
	int numData = 10;

	if (axis == true)	numData++;
	if (time_integration != 0) numData++;
	if (viscous == true) numData += 2;

	Data::DataStorageVector<Data::DataStorage>	sampleData(numData);

	numData = 0;
	sampleData(numData)	= CFD_VariableSet::Q(species_set, ND, NER, NEV, NEE, viscous, axis);		numData++;
	sampleData(numData)	= CFD_VariableSet::V(species_set, ND, NER, NEV, NEE, viscous, axis);		numData++;
	sampleData(numData)	= CFD_VariableSet::W(species_set, ND, NER, NEV, NEE, viscous, axis);		numData++;
	sampleData(numData)	= CFD_VariableSet::MIX(species_set, ND, NER, NEV, NEE, viscous, axis);		numData++;
	sampleData(numData)	= CFD_VariableSet::Xs(species_set, ND, NER, NEV, NEE, viscous, axis);		numData++;
	sampleData(numData)	= CFD_VariableSet::Ys(species_set, ND, NER, NEV, NEE, viscous, axis);		numData++;
	sampleData(numData)	= CFD_VariableSet::Residue(species_set, ND, NER, NEV, NEE, viscous, axis);	numData++;
	sampleData(numData)	= CFD_VariableSet::dQ(species_set, ND, NER, NEV, NEE, viscous, axis); 		numData++;
	sampleData(numData)	= CFD_VariableSet::Qnew(species_set, ND, NER, NEV, NEE, viscous, axis); 	numData++;
	sampleData(numData)	= CFD_VariableSet::Source(species_set, ND, NER, NEV, NEE, viscous, axis);	numData++;

	if (axis == true)
	{
		Common::Map1D<std::string, int>	divVcMap (1);
		divVcMap.insert("divergence V", 0);
		Data::DataStorage data_CFD_divVc(NAME_DIVV, 1, divVcMap);

		sampleData(numData)	= data_CFD_divVc; numData++;
	}

	if (time_integration != 0)
	{
		sampleData(numData)	= CFD_VariableSet::dpdQ(species_set, ND, NER, NEV, NEE, viscous, axis);	numData++;
	}

	if (viscous == true)
	{
		sampleData(numData)	= CFD_VariableSet::DiffusionCoeff(species_set, ND, NER, NEV, NEE, viscous, axis);	numData++;
		sampleData(numData)	= CFD_VariableSet::ViscosityCoeff(species_set, ND, NER, NEV, NEE, viscous, axis);	numData++;

	}

	sampleData.mapping();


	return(sampleData);
}






Data::DataStorageVector<Data::DataStorage2D> CFD_DataTemplateWIG::CellData2D(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis, int time_integration)
{
	int numData = 0;

	if (time_integration != 0) numData++;
	if (viscous == true) numData ++;


	Data::DataStorageVector<Data::DataStorage2D>	sampleData(numData);

	numData = 0;
	if (time_integration != 0)
	{
		sampleData(numData)	= CFD_VariableSet::dTdQ(species_set, ND, NER, NEV, NEE, viscous, axis);	numData++;
	}

	if (viscous == true)
	{
		sampleData(numData)	= CFD_VariableSet::thermal_conductivity(species_set, ND, NER, NEV, NEE, viscous, axis);	numData++;
	}

	sampleData.mapping();


	return(sampleData);
}








}
}
