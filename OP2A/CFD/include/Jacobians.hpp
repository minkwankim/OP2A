/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 10, 2015
 *      			Author: Minkwan Kim
 *
 * Jacobians.hpp
 * 			-  
 *  
 */



#include <limits>
#include "Common/include/Exception_NaNValue.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"

#include "CFD/include/VariableChange.hpp"
#include "CFD/include/VariableConstants.hpp"
#include "CFD/include/Derivatives.hpp"


namespace OP2A{
namespace CFD{


class CFD_API Jacobians: public Common::NonInstantiable<Jacobians>
{
public:
	static void Axisymmertic(Data::DataStorage& dpdQ, int NS, int ND, int NE, double S, Data::DataStorage2D& dSdQ);
	static void S_he(Data::DataStorageVector<Data::DataStorage>& data1D, Data::DataStorage2D& dT, CHEM::SpeciesSet& species_set, int ND,
			 	 	 unsigned int indexTe, unsigned int indexQ, unsigned int indexV, unsigned int indexWint,
			 	 	 double S, Data::DataStorage& dS);


};


}
}
