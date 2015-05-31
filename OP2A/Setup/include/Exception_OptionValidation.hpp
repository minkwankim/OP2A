/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 26, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_OptionValidation.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_OPTIONVALIDATION_HPP_
#define EXCEPTION_OPTIONVALIDATION_HPP_


#include	"Common/include/Exception.hpp"

#include	"Setup/include/SetupAPI.hpp"

namespace	OP2A{
namespace	Setup{

class Setup_API	ExceptionOptionValidation	: public Common::Exception
{
public:
	ExceptionOptionValidation(const Common::Code_location& where, const std::string& what) : Common::Exception(where, what,	"ExceptionOptionValidation")
	{

	}


	ExceptionOptionValidation(const ExceptionOptionValidation& e) throw() : Common::Exception(e)
	{

	}
};


}
}


#endif /* EXCEPTION_OPTIONVALIDATION_HPP_ */
