/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 21, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_BadMatch.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_BADMATCH_HPP_
#define EXCEPTION_BADMATCH_HPP_

#include	"Common/include/Exception.hpp"

#include	"Setup/include/SetupAPI.hpp"


/*
 * Class for exception handling
 * @author	Minkwan Kim
 * @version	1.0	20/05/2015
 */
namespace	OP2A{
namespace	Setup{

class Setup_API	ExceptionBadMatch	: public Common::Exception
{
public:
	ExceptionBadMatch(const Common::Code_location& where, const std::string& what) : Common::Exception(where, what,	"ExceptionBasMatch")
	{

	}


	ExceptionBadMatch(const ExceptionBadMatch& e) throw() : Common::Exception(e)
	{

	}
};


}
}


#endif /* EXCEPTION_BADMATCH_HPP_ */
