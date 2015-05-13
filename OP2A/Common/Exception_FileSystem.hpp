/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_FileSystem.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_FILESYSTEM_HPP_
#define EXCEPTION_FILESYSTEM_HPP_

#include "Common/include/Exception.hpp"

namespace OP2A{
namespace Common{

class Common_API ExceptionFileSystem : public Common::Exception
{

	ExceptionFileSystem ( const Common::Code_location& where, const std::string& what)
			: Common::Exception(where,what,"ExceptionFileSystem")
	{

	}


	ExceptionFileSystem (const ExceptionFileSystem& e) throw () : Common::Exception(e)
	{

	}
};


}
}



#endif /* EXCEPTION_FILESYSTEM_HPP_ */
