/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 18, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChangeCommon.hpp
 * 			-  
 *  
 */
#ifndef VARIABLECHANGECOMMON_HPP_
#define VARIABLECHANGECOMMON_HPP_

#include "CFD/include/CFDCommon.hpp"


namespace OP2A{
namespace CFD{



/*
 * Class for Common module for variable change
 * 	- variable change related to mass and momentum parts
 * 	@author Minkwan
 * 	@version 1.0 18/6/2015
 */

class CFD_API VariableChangeCommon : public Common::NonInstantiable<VariableChangeCommon>
{
public:
	static void V_to_Q(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT);
	static void Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT);

	static void V_to_W(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT);
	static void W_to_V(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT);

	static void Q_to_W(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT);
	static void W_to_Q(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT);
};



}
}



#endif /* VARIABLECHANGECOMMON_HPP_ */
