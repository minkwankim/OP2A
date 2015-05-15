/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_FloatingPoint.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_FLOATINGPOINT_HPP_
#define EXCEPTION_FLOATINGPOINT_HPP_

#include "Common/include/Exception.hpp"

namespace OP2A{
namespace Common{

class Common_API ExceptionFloatingPoint: public Common::Exception
{
public:
	ExceptionFloatingPoint ( const Common::Code_location& where, const std::string& what)
			: Common::Exception(where,what,"ExceptionFloatingPoint")
	{

	}


	ExceptionFloatingPoint (const ExceptionFloatingPoint& e) throw () : Common::Exception(e)
	{

	}
};


}
}



#endif /* EXCEPTION_FLOATINGPOINT_HPP_ */
