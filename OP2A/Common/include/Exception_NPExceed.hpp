/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 25, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_NPExceed.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_NPEXCEED_HPP_
#define EXCEPTION_NPEXCEED_HPP_


#include "Common/include/Exception.hpp"

namespace OP2A{
namespace Common{

class Common_API ExceptionNPExceed: public Common::Exception
{
public:

	ExceptionNPExceed ( const Common::Code_location& where, const std::string& what)
			: Common::Exception(where, what, "ExceptionNPExceed")
	{

	}


	ExceptionNPExceed (const ExceptionNPExceed& e) throw () : Common::Exception(e)
	{

	}
};


}
}



#endif /* EXCEPTION_NPEXCEED_HPP_ */
