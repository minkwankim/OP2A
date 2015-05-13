/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_NoSuchValue.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_NOSUCHVALUE_HPP_
#define EXCEPTION_NOSUCHVALUE_HPP_

#include "Common/include/Exception.hpp"

namespace OP2A{
namespace Common{

class Common_API ExceptionNoSuchValue: public Common::Exception
{
public:

	ExceptionNoSuchValue ( const Common::Code_location& where, const std::string& what)
			: Common::Exception(where,what,"ExceptionNoSuchValue")
	{

	}


	ExceptionNoSuchValue (const ExceptionNoSuchValue& e) throw () : Common::Exception(e)
	{

	}
};


}
}

#endif /* EXCEPTION_NOSUCHVALUE_HPP_ */
