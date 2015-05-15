/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_NoSuchStorage.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_NOSUCHSTORAGE_HPP_
#define EXCEPTION_NOSUCHSTORAGE_HPP_



#include "Common/include/Exception.hpp"

namespace OP2A{
namespace Common{

class Common_API ExceptionNoSuchStorage: public Common::Exception
{
public:

	ExceptionNoSuchStorage ( const Common::Code_location& where, const std::string& what)
			: Common::Exception(where,what,"ExceptionNoSuchStorage")
	{

	}


	ExceptionNoSuchStorage (const ExceptionNoSuchStorage& e) throw () : Common::Exception(e)
	{

	}
};


}
}


#endif /* EXCEPTION_NOSUCHSTORAGE_HPP_ */
