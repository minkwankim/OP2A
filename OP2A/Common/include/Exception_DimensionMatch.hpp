/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 5, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_DimensionMatch.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_DIMENSIONMATCH_HPP_
#define EXCEPTION_DIMENSIONMATCH_HPP_



#include "Common/include/Exception.hpp"

namespace OP2A {

namespace Common {

class Common_API ExceptionDimensionMatch : public Common::Exception
{
public:
	ExceptionDimensionMatch (const Common::Code_location& where, const std::string& what) :
    Common::Exception(where, what, "ExceptionDimensionMatch") {}

  /// Copy constructor
	ExceptionDimensionMatch(const ExceptionDimensionMatch& e) throw  () : Exception(e) {}
};


}
}

#endif /* EXCEPTION_DIMENSIONMATCH_HPP_ */
