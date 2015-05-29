/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 3, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_DataStorageVector.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_DATASTORAGEVECTOR_HPP_
#define EXCEPTION_DATASTORAGEVECTOR_HPP_



#include "Common/include/Exception.hpp"

namespace OP2A{
namespace Data{

class ExceptionDataStorageVector :public Common::Exception
{
public:
	ExceptionDataStorageVector (const Common::Code_location& where, const std::string& what) :
		Common::Exception(where,what,"ExceptionDataStorageVector")
	{

	}

	ExceptionDataStorageVector(const ExceptionDataStorageVector& e) throw() : Common::Exception(e)
	{
	}
};


}
}



#endif /* EXCEPTION_DATASTORAGEVECTOR_HPP_ */
