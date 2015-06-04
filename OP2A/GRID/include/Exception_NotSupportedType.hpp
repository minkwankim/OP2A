/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 4, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_NotSupportedType.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_NOTSUPPORTEDTYPE_HPP_
#define EXCEPTION_NOTSUPPORTEDTYPE_HPP_


#include "Common/include/Exception.hpp"

namespace OP2A{
namespace GRID{

class ExceptionNotSupportedType :public Common::Exception
{
public:
	ExceptionNotSupportedType (const Common::Code_location& where, const std::string& what) :
		Common::Exception(where,what,"ExceptionNotSupportedType")
	{

	}

	ExceptionNotSupportedType(const ExceptionNotSupportedType& e) throw() : Common::Exception(e)
	{
	}
};


}
}



#endif /* EXCEPTION_NOTSUPPORTEDTYPE_HPP_ */
