/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_LibLoader.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_LIBLOADER_HPP_
#define EXCEPTION_LIBLOADER_HPP_



#include "Common/include/Exception.hpp"

namespace OP2A{
namespace Common{

class ExceptionLibLoader :public Common::Exception
{
public:
	ExceptionLibLoader (const Common::Code_location& where, const std::string& what) :
		Common::Exception(where,what,"ExceptionLibLoader")
	{

	}

	ExceptionLibLoader(const ExceptionLibLoader& e) throw() : Common::Exception(e)
	{
	}
};


}
}

#endif /* EXCEPTION_LIBLOADER_HPP_ */
