/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * DemangledTypeID.hpp
 * 			-  
 *  
 */
#ifndef DEMANGLEDTYPEID_HPP_
#define DEMANGLEDTYPEID_HPP_

#include "Common/include/Common.hpp"
#include "Common/include/CommonAPI.hpp"

namespace OP2A{
namespace Common{

Common_API std::string demangle(const char* type);
#define DEMANGLED_TYPEID(a) OP2A::Common::demangle(typeid(a).name())




}
}


#endif /* DEMANGLEDTYPEID_HPP_ */
