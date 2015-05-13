/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_FailedCast.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_FAILEDCAST_HPP_
#define EXCEPTION_FAILEDCAST_HPP_


#include "Common/include/Exception.hpp"


namespace OP2A{
namespace Common{

class ExceptionFailedCast : public Common::Exception
{
public:
	ExceptionFailedCast ( const Common::Code_location& where, const std::string& what)
			:Common::Exception(where, what, "ExceptionFailedCast") {}

	ExceptionFailedCast ( const ExceptionFailedCast& e) throw () : Exception(e) {}
};

}
}

#endif /* EXCEPTION_FAILEDCAST_HPP_ */
