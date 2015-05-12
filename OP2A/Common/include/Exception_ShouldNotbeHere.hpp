/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_ShouldNotbeHere.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_SHOULDNOTBEHERE_HPP_
#define EXCEPTION_SHOULDNOTBEHERE_HPP_


#include "Common/include/Exception.hpp"


namespace OP2A{
namespace Common{

class Common_API ExceptionShouldNotbeHere: public Common::Exception
{

	ExceptionShouldNotbeHere ( const Common::Code_location& where, const std::string& what)
			: Common::Exception(where, what, "ExceptionShouldNotbeHere")
	{

	}


	ExceptionShouldNotbeHere (const ExceptionShouldNotbeHere& e) throw () : Common::Exception(e)
	{

	}
};


}
}



#endif /* EXCEPTION_SHOULDNOTBEHERE_HPP_ */
