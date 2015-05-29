/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 29, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_DataStorageSize.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_DATASTORAGESIZE_HPP_
#define EXCEPTION_DATASTORAGESIZE_HPP_




#include "Common/include/Exception.hpp"

namespace OP2A{
namespace Data{

class ExceptionDataStorageSize :public Common::Exception
{
public:
	ExceptionDataStorageSize (const Common::Code_location& where, const std::string& what) :
		Common::Exception(where,what,"ExceptionDataStorageSize")
	{

	}

	ExceptionDataStorageSize(const ExceptionDataStorageSize& e) throw() : Common::Exception(e)
	{
	}
};


}
}


#endif /* EXCEPTION_DATASTORAGESIZE_HPP_ */
