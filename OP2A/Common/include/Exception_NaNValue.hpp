/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 5, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_NANValue.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_NANVALUE_HPP_
#define EXCEPTION_NANVALUE_HPP_


namespace OP2A{
namespace Common{

class ExceptionNaNValue :public Common::Exception
{
public:
	ExceptionNaNValue (const Common::Code_location& where, const std::string& what) :
		Common::Exception(where,what,"ExceptionNaNValue")
	{

	}

	ExceptionNaNValue(const ExceptionNaNValue& e) throw() : Common::Exception(e)
	{
	}
};


}
}




#endif /* EXCEPTION_NANVALUE_HPP_ */
