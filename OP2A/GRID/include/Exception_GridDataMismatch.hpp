/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 4, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_GridDataMismatch.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_GRIDDATAMISMATCH_HPP_
#define EXCEPTION_GRIDDATAMISMATCH_HPP_



#include "Common/include/Exception.hpp"

namespace OP2A{
namespace GRID{

class ExceptionGridDataMismatch :public Common::Exception
{
public:
	ExceptionGridDataMismatch (const Common::Code_location& where, const std::string& what) :
		Common::Exception(where,what,"ExceptionGridDataMismatch")
	{

	}

	ExceptionGridDataMismatch(const ExceptionGridDataMismatch& e) throw() : Common::Exception(e)
	{
	}
};


}
}



#endif /* EXCEPTION_GRIDDATAMISMATCH_HPP_ */
