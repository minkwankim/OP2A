/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 25, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_OutOfRange.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_OUTOFRANGE_HPP_
#define EXCEPTION_OUTOFRANGE_HPP_

#include "Common/include/Exception.hpp"

namespace OP2A{
namespace Common{

class Common_API ExceptionOutOfRange: public Common::Exception
{
public:

	ExceptionOutOfRange ( const Common::Code_location& where, const std::string& what)
			: Common::Exception(where, what, "ExceptionOutOfRange")
	{

	}


	ExceptionOutOfRange (const ExceptionOutOfRange& e) throw () : Common::Exception(e)
	{

	}
};


}
}



#endif /* EXCEPTION_OUTOFRANGE_HPP_ */
