/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 17, 2015
 *      			Author: Minkwan Kim
 *
 * VariableSamples.hpp
 * 			-  
 *  
 */
#ifndef VARIABLESAMPLES_HPP_
#define VARIABLESAMPLES_HPP_


#include "CFD/include/VariableSet.hpp"


namespace OP2A{
namespace CFD{

class CFD_API 	CFD_DataTemplateWIG: public Common::NonInstantiable<CFD_VariableSet>
{
public:
	static Data::DataStorageVector<Data::DataStorage>		CellData1D(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis, int time_integration);
	static Data::DataStorageVector<Data::DataStorage2D>		CellData2D(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis, int time_integration);

	static Data::DataStorageVector<Data::DataStorage>		FaceData1D(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis, int time_integration);
	static Data::DataStorageVector<Data::DataStorage2D>		FaceData2D(const CHEM::SpeciesSet& species_set, int ND, int NER, int NEV, int NEE, bool viscous, bool axis, int time_integration);
};



}
}

#endif /* VARIABLESAMPLES_HPP_ */
