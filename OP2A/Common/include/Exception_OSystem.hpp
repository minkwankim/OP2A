/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_OSystem.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_OSYSTEM_HPP_
#define EXCEPTION_OSYSTEM_HPP_



#include "Common/include/Exception.hpp"

namespace OP2A{
namespace Common{

class Common_API ExceptionOSystem: public Common::Exception
{
public:

	ExceptionOSystem ( const Common::Code_location& where, const std::string& what)
			: Common::Exception(where, what, "ExceptionOSystem")
	{

	}


	ExceptionOSystem (const ExceptionOSystem& e) throw () : Common::Exception(e)
	{

	}
};


}
}

#endif /* EXCEPTION_OSYSTEM_HPP_ */
