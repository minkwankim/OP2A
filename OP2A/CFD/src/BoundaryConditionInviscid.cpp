/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 23, 2015
 *      			Author: Minkwan Kim
 *
 * BoundaryConditionInviscid.hpp
 * 			-  
 *  
 */




#include "CFD/include/BoundaryConditions.hpp"

namespace OP2A{
namespace CFD{


CFDBCTypes	BCInviscid::BCTypeInCFD(const int FaceBCType)
{
	CFDBCTypes	bctype;

	switch (FaceBCType)
	{
	case GRID::BCType::interior:
		bctype	= CFDBCTypes::Interior;
		break;

	case GRID::BCType::wall:
		bctype	= CFDBCTypes::WallType;
		break;

	case GRID::BCType::inlet:
		bctype	= CFDBCTypes::InletType;
		break;

	case GRID::BCType::outlet:
		bctype	= CFDBCTypes::ExitType;
		break;

	case GRID::BCType::freestream:
		bctype	= CFDBCTypes::FreestreamType;
		break;

	case GRID::BCType::symmetric:
		bctype	= CFDBCTypes::WallType;
		break;

	case GRID::BCType::axis:
		bctype	= CFDBCTypes::WallType;
		break;

	case GRID::BCType::anode:
		bctype	= CFDBCTypes::WallType;
		break;

	case GRID::BCType::cathode:
		bctype	= CFDBCTypes::WallType;
		break;

	case GRID::BCType::dielectricwall:
		bctype	= CFDBCTypes::WallType;
		break;
	}

	return bctype;
}







}
}
