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

#include "CFD/include/OP2A_CFD.hpp"


namespace OP2A{
namespace CFD{



/*
 * NER = 0
 * NEV = 0
 * NEE = 0
 */

class CFD_API VariableChangeCommon : public Common::NonInstantiable<VariableChangeCommon>
{
public:
	static void V_to_Q(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q);
	static void Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V);

	static void V_to_W(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W);
	static void W_to_V(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V);

	static void Q_to_W(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W);
	static void W_to_Q(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q);


};



}
}



#endif /* VARIABLECHANGECOMMON_HPP_ */
