/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 18, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChange.hpp
 * 			-  
 *  
 */
#ifndef VARIABLECHANGE_HPP_
#define VARIABLECHANGE_HPP_

#include "VariableSet.hpp"
#include "VariableConstants.hpp"

#include "CFD/include/VariableChangeCommon.hpp"
#include "CFD/include/VariableChangeType1.hpp"


namespace OP2A{
namespace CFD{


class CFD_API VariableChange : public Common::NonInstantiable<VariableChange>
{
public:
	static void V_to_Q(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q);
	static void Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V);

	static void V_to_W(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W);
	static void W_to_V(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V);

	static void Q_to_W(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W);
	static void W_to_Q(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V);
};


}
}

#endif /* VARIABLECHANGE_HPP_ */
