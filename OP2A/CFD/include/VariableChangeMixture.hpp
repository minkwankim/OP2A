/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 16, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChangeMixture.hpp
 * 			-  
 *  
 */


#ifndef VARIABLECHANGEMIXTURE_HPP_
#define VARIABLECHANGEMIXTURE_HPP_

#include <limits>
#include "CFD/include/OP2A_CFD.hpp"


namespace OP2A{
namespace CFD{


class CFD_API VariableChangeMixture : public Common::NonInstantiable<VariableChangeMixture>
{
public:
	void Xs(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Xs);
	void Ys(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Xs);
};


}
}



#endif /* VARIABLECHANGEMIXTURE_HPP_ */
