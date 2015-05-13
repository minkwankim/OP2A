/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_FileFormat.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_FILEFORMAT_HPP_
#define EXCEPTION_FILEFORMAT_HPP_


#include "Common/include/Exception.hpp"

namespace OP2A{
namespace Common{

class ExceptionFileFormat :public Common::Exception
{
public:
	ExceptionFileFormat (const Common::Code_location& where, const std::string& what) :
		Common::Exception(where,what,"ExceptionFileFormat")
	{

	}

	ExceptionFileFormat(const ExceptionFileFormat& e) throw() : Common::Exception(e)
	{
	}
};


}
}


#endif /* EXCEPTION_FILEFORMAT_HPP_ */
