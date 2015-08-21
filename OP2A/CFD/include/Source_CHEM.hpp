/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 17, 2015
 *      			Author: Minkwan Kim
 *
 * Source_CHEM.hpp
 * 			-  
 *  
 */
#ifndef SOURCE_CHEM_HPP_
#define SOURCE_CHEM_HPP_



#include "CFD/include/CFDCommon.hpp"
#include "DATA/include/DataStorage2D.hpp"
#include "DATA/include/DataStorageVector.hpp"



namespace OP2A{
namespace CFD{


class CFD_API Source_CHEM : public Common::NonInstantiable<Source_CHEM>
{
public:
	static void S_chem(Data::DataStorageVector<Data::DataStorage>& data1D, CHEM::SpeciesSet& species_set, int ND,
					   int type, unsigned int indexQ, unsigned int indexV, unsigned int indexW,
					   Data::DataStorage& Sc);
};





}
}


#endif /* SOURCE_CHEM_HPP_ */
