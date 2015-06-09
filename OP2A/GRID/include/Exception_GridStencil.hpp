/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 9, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_GridStencil.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_GRIDSTENCIL_HPP_
#define EXCEPTION_GRIDSTENCIL_HPP_



#include "Common/include/Exception.hpp"

namespace OP2A{
namespace GRID{

class ExceptionGridStencil :public Common::Exception
{
public:
	ExceptionGridStencil (const Common::Code_location& where, const std::string& what) :
		Common::Exception(where,what,"ExceptionGridStencil")
	{

	}

	ExceptionGridStencil(const ExceptionGridStencil& e) throw() : Common::Exception(e)
	{
	}
};


}
}

#endif /* EXCEPTION_GRIDSTENCIL_HPP_ */
