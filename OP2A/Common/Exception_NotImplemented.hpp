/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_NotImplemented.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_NOTIMPLEMENTED_HPP_
#define EXCEPTION_NOTIMPLEMENTED_HPP_


#include "Common/include/Exception.hpp"

namespace OP2A{
namespace Common{

class Common_API ExceptionNotImplemented: public Common::Exception
{

	ExceptionNotImplemented ( const Common::Code_location& where, const std::string& what)
			: Common::Exception(where, what, "ExceptionNotImplemented")
	{

	}


	ExceptionNotImplemented (const ExceptionNotImplemented& e) throw () : Common::Exception(e)
	{

	}
};


}
}



#endif /* EXCEPTION_NOTIMPLEMENTED_HPP_ */
