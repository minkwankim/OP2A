/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_Parallel.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_PARALLEL_HPP_
#define EXCEPTION_PARALLEL_HPP_

#include "Common/include/Exception.hpp"

namespace OP2A{
namespace Common{

class Common_API ExceptionParallel: public Common::Exception
{

	ExceptionParallel ( const Common::Code_location& where, const std::string& what)
			: Common::Exception(where, what, "ExceptionParallel")
	{

	}


	ExceptionParallel (const ExceptionParallel& e) throw () : Common::Exception(e)
	{

	}
};


}
}




#endif /* EXCEPTION_PARALLEL_HPP_ */
