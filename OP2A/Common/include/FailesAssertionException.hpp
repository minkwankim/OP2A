/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 11, 2015
 *      			Author: Minkwan Kim
 *
 * FailesAssertionException.hpp
 * 			-  
 *  
 */
#ifndef FAILESASSERTIONEXCEPTION_HPP_
#define FAILESASSERTIONEXCEPTION_HPP_


#include "Common/include/Exception.hpp"

//////////////////////////////////////////////////////////////////////////////

namespace OP2A {

namespace Common {

class Common_API FailedAssertionException : public Common::Exception {
public:

  FailedAssertionException (const Common::Code_location& where, const std::string& what) :
    Common::Exception(where, what, "FailedAssertionException") {}

  /// Copy constructor
  FailedAssertionException(const FailedAssertionException& e) throw  () : Exception(e) {}

};


}
}
#endif /* FAILESASSERTIONEXCEPTION_HPP_ */
