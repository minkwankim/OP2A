/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 15, 2015
 *      			Author: Minkwan Kim
 *
 * Source_NONEQ.hpp
 * 			-  
 *  
 */
#ifndef SOURCE_NONEQ_HPP_
#define SOURCE_NONEQ_HPP_


#include "CFD/include/CFDCommon.hpp"
#include "DATA/include/DataStorage2D.hpp"
#include "DATA/include/DataStorageVector.hpp"



namespace OP2A{
namespace CFD{


class CFD_API Source_NONEQ : public Common::NonInstantiable<Source_NONEQ>
{
public:
	static double S_he(Data::DataStorageVector<Data::DataStorage>& data1D, CHEM::SpeciesSet& species_set, int ND,
					   unsigned int indexTe, unsigned int indexQ, unsigned int indexV, unsigned int indexW);


	static void S_NONEQ(Data::DataStorageVector<Data::DataStorage>& data1D, CHEM::SpeciesSet& species_set, int ND,
						unsigned int type, unsigned int indexQ, unsigned int indexV, unsigned int indexW, unsigned int indexS);

};





}
}


#endif /* SOURCE_NONEQ_HPP_ */
