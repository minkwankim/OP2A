/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 16, 2015
 *      			Author: Minkwan Kim
 *
 * ErrorCheckNonNegative.hpp
 * 			-  
 *  
 */
#ifndef ERRORCHECKNONNEGATIVE_HPP_
#define ERRORCHECKNONNEGATIVE_HPP_

#include <limits>
#include "Common/include/Common.hpp"
#include "Common/include/Exception_NaNValue.hpp"
#include "Common/include/Exception_NegativeValue.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"


namespace OP2A {
namespace Common {


template <typename T>
class ErrorCheckNonNegative
{
public:
	explicit ErrorCheckNonNegative(T value, const std::string& moduleName)
	{
		std::ostringstream message;

		if (value	!= value)
		{
			message << "Nan value for " << moduleName;
			throw Common::ExceptionNaNValue (FromHere(), message.str());
		}

		if (value == std::numeric_limits<double>::infinity())
		{
			message << "Inifnite Value for " << moduleName;
			throw Common::ExceptionInfiniteValue (FromHere(), message.str());
		}

		if (value < 0)
		{
			message << "Negative Value for " << moduleName;
			throw Common::ExceptionNegativeValue (FromHere(), message.str());
		}
	}
};



}
}
#endif /* ERRORCHECKNONNEGATIVE_HPP_ */
