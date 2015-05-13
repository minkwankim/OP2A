/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_Parser.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_PARSER_HPP_
#define EXCEPTION_PARSER_HPP_

#include "Common/include/Exception.hpp"


namespace OP2A{
namespace Common{

class Common_API ExceptionParser: public Common::Exception
{

	ExceptionParser ( const Common::Code_location& where, const std::string& what)
			: Common::Exception(where, what, "ExceptionParser")
	{

	}


	ExceptionParser (const ExceptionParser& e) throw () : Common::Exception(e)
	{

	}
};


}
}

#endif /* EXCEPTION_PARSER_HPP_ */
