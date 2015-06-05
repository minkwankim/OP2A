/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 5, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_InfiniteValue.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_INFINITEVALUE_HPP_
#define EXCEPTION_INFINITEVALUE_HPP_

namespace OP2A{
namespace Common{

class ExceptionInfiniteValue :public Common::Exception
{
public:
	ExceptionInfiniteValue (const Common::Code_location& where, const std::string& what) :
		Common::Exception(where,what,"ExceptionInfiniteValue")
	{

	}

	ExceptionInfiniteValue(const ExceptionInfiniteValue& e) throw() : Common::Exception(e)
	{
	}
};


}
}



#endif /* EXCEPTION_INFINITEVALUE_HPP_ */
