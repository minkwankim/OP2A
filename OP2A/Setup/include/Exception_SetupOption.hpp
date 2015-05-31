/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 26, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_SetupOption.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_SETUPOPTION_HPP_
#define EXCEPTION_SETUPOPTION_HPP_

#include	"Common/include/Exception.hpp"

#include	"Setup/include/SetupAPI.hpp"

namespace	OP2A{
namespace	Setup{

class Setup_API	ExceptionSetupOption	: public Common::Exception
{
public:
	ExceptionSetupOption(const Common::Code_location& where, const std::string& what) : Common::Exception(where, what,	"ExceptionSetupOption")
	{

	}


	ExceptionSetupOption(const ExceptionSetupOption& e) throw() : Common::Exception(e)
	{

	}
};


}
}




#endif /* EXCEPTION_SETUPOPTION_HPP_ */
