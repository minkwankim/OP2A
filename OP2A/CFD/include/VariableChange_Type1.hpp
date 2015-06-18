/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 15, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChange_V_to_Q.hpp
 * 			-  
 *  
 */
#ifndef VARIABLECHANGE_V_TO_Q_HPP_
#define VARIABLECHANGE_V_TO_Q_HPP_


#include "CFD/include/OP2A_CFD.hpp"


namespace OP2A{
namespace CFD{



/*
 * NER = 0
 * NEV = 0
 * NEE = 0
 */

class CFD_API VariableChangeType1 : public Common::NonInstantiable<VariableChangeType1>
{
public:
	void V_to_Q(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q);
	void Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V);


};



}
}


#endif /* VARIABLECHANGE_V_TO_Q_HPP_ */
