/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_Setup.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_SETUP_HPP_
#define EXCEPTION_SETUP_HPP_

#include "Common/include/Exception.hpp"

namespace OP2A{
namespace Common{

class Common_API ExceptionSetup: public Common::Exception
{
public:

	ExceptionSetup ( const Common::Code_location& where, const std::string& what)
			: Common::Exception(where, what, "ExceptionSetup")
	{

	}


	ExceptionSetup (const ExceptionSetup& e) throw () : Common::Exception(e)
	{

	}
};


}
}


#endif /* EXCEPTION_SETUP_HPP_ */
