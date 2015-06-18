/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 16, 2015
 *      			Author: Minkwan Kim
 *
 * OPPA_CFD.hpp
 * 			-  
 *  
 */
#ifndef OPPA_CFD_HPP_
#define OPPA_CFD_HPP_

#include <limits>

#include "Common/include/ErrorCheckNonNegative.hpp"

#include "DATA/include/DataStorage.hpp"
#include "DATA/include/DataStorage2D.hpp"
#include "DATA/include/DataStorageVector.hpp"

#include "CHEM/include/SpeciesSet.hpp"

#include "CFD/include/CFD_API.hpp"


namespace OP2A{
namespace CFD{

static unsigned int CFD_type;
static unsigned int CFD_NT;


}
}


#endif /* OPPA_CFD_HPP_ */
