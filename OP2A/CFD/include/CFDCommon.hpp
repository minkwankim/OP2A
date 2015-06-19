/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 18, 2015
 *      			Author: Minkwan Kim
 *
 * CFDCommon.hpp
 * 			-  
 *  
 */
#ifndef CFDCOMMON_HPP_
#define CFDCOMMON_HPP_


#include <limits>
#include "Common/include/ErrorCheckNonNegative.hpp"

#include "DATA/include/DataStorage.hpp"
#include "DATA/include/DataStorage2D.hpp"
#include "DATA/include/DataStorageVector.hpp"

#include "CHEM/include/SpeciesSet.hpp"


#include "CFD/include/CFD_API.hpp"
#include "CFD/include/VariableConstants.hpp"

namespace OP2A{
namespace CFD{

double CFD_calculate_Tve(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, double Eve, unsigned int iterMax, double eps0, unsigned int CFD_NT);


}
}

#endif /* CFDCOMMON_HPP_ */
