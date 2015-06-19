/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 19, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_CatchInLoop.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_CATCHINLOOP_HPP_
#define EXCEPTION_CATCHINLOOP_HPP_


#include "Common/include/Exception.hpp"

namespace OP2A {

namespace Common {

class Common_API ExceptionCatchInLoop : public Common::Exception
{
public:
	ExceptionCatchInLoop (const Common::Code_location& where, const std::string& what) :
    Common::Exception(where, what, "ExceptionCatchInLoop") {}

  /// Copy constructor
	ExceptionCatchInLoop(const ExceptionCatchInLoop& e) throw  () : Exception(e) {}
};


}
}


#endif /* EXCEPTION_CATCHINLOOP_HPP_ */
