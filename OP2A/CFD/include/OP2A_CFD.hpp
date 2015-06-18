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


#include "CFD/include/VariableSet.hpp"


namespace OP2A{
namespace CFD{

static unsigned int CFD_NT   = 1;
static unsigned variableType = 0;

void assignVariableType(int NER, int NEV, int NEE);



}
}




#endif /* OPPA_CFD_HPP_ */
