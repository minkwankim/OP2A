/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 14, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_BadValue.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_BADVALUE_HPP_
#define EXCEPTION_BADVALUE_HPP_


#include "Common/include/Exception.hpp"

//////////////////////////////////////////////////////////////////////////////

namespace OP2A {

namespace Common {

class Common_API ExceptionBadValue : public Common::Exception
{
public:
	ExceptionBadValue (const Common::Code_location& where, const std::string& what) :
    Common::Exception(where, what, "ExceptionBadValue") {}

	/// Copy constructor
	ExceptionBadValue(const ExceptionBadValue& e) throw  () : Exception(e) {}
};


}
}


#endif /* EXCEPTION_BADVALUE_HPP_ */
