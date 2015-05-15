/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 14, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_MempryAllocator.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_MEMPRYALLOCATOR_HPP_
#define EXCEPTION_MEMPRYALLOCATOR_HPP_


#include "Common/include/Exception.hpp"

//////////////////////////////////////////////////////////////////////////////

namespace OP2A {

namespace Common {

class Common_API ExceptionMemoryAllocator : public Common::Exception
{
public:
	ExceptionMemoryAllocator (const Common::Code_location& where, const std::string& what) :
    Common::Exception(where, what, "ExceptionMemoryAllocator") {}

	/// Copy constructor
	ExceptionMemoryAllocator(const ExceptionMemoryAllocator& e) throw  () : Exception(e) {}
};


}
}


#endif /* EXCEPTION_MEMPRYALLOCATOR_HPP_ */
