/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 4, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_DimensionMisMatch.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_DIMENSIONMISMATCH_HPP_
#define EXCEPTION_DIMENSIONMISMATCH_HPP_


#include "Common/include/Exception.hpp"

namespace OP2A{
namespace GRID{

class ExceptionDimensionMismatch :public Common::Exception
{
public:
	ExceptionDimensionMismatch (const Common::Code_location& where, const std::string& what) :
		Common::Exception(where,what,"ExceptionDimensionMismatch")
	{

	}

	ExceptionDimensionMismatch(const ExceptionDimensionMismatch& e) throw() : Common::Exception(e)
	{
	}
};


}
}


#endif /* EXCEPTION_DIMENSIONMISMATCH_HPP_ */
