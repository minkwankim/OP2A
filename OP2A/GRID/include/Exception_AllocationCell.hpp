/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 4, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_AllocationCell.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_ALLOCATIONCELL_HPP_
#define EXCEPTION_ALLOCATIONCELL_HPP_


#include "Common/include/Exception.hpp"

namespace OP2A{
namespace GRID{

class ExceptionAllocationCell :public Common::Exception
{
public:
	ExceptionAllocationCell (const Common::Code_location& where, const std::string& what) :
		Common::Exception(where,what,"ExceptionAllocationCell")
	{

	}

	ExceptionAllocationCell(const ExceptionAllocationCell& e) throw() : Common::Exception(e)
	{
	}
};


}
}




#endif /* EXCEPTION_ALLOCATIONCELL_HPP_ */
