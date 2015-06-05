/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 5, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_CellConstruction.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_CELLCONSTRUCTION_HPP_
#define EXCEPTION_CELLCONSTRUCTION_HPP_



namespace OP2A{
namespace GRID{

class ExceptionCellConstruction :public Common::Exception
{
public:
	ExceptionCellConstruction (const Common::Code_location& where, const std::string& what) :
		Common::Exception(where,what,"ExceptionCellConstruction")
	{

	}

	ExceptionCellConstruction(const ExceptionCellConstruction& e) throw() : Common::Exception(e)
	{
	}
};


}
}

#endif /* EXCEPTION_CELLCONSTRUCTION_HPP_ */
