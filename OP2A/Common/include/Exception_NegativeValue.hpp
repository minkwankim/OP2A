/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 5, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_NagativeValue.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_NAGATIVEVALUE_HPP_
#define EXCEPTION_NAGATIVEVALUE_HPP_


namespace OP2A{
namespace Common{

class ExceptionNegativeValue :public Common::Exception
{
public:
	ExceptionNegativeValue (const Common::Code_location& where, const std::string& what) :
		Common::Exception(where,what,"ExceptionNegativeValue")
	{

	}

	ExceptionNegativeValue(const ExceptionNegativeValue& e) throw() : Common::Exception(e)
	{
	}
};


}
}



#endif /* EXCEPTION_NAGATIVEVALUE_HPP_ */
