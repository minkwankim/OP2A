/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_NullPointer.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_NULLPOINTER_HPP_
#define EXCEPTION_NULLPOINTER_HPP_

#include "Common/include/Exception.hpp"

namespace OP2A{
namespace Common{

class Common_API ExceptionNullPointer: public Common::Exception
{
public:

	ExceptionNullPointer ( const Common::Code_location& where, const std::string& what)
			: Common::Exception(where, what, "ExceptionNullPointer")
	{

	}


	ExceptionNullPointer (const ExceptionNullPointer& e) throw () : Common::Exception(e)
	{

	}
};


}
}


#endif /* EXCEPTION_NULLPOINTER_HPP_ */
