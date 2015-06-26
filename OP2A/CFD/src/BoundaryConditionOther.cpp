/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 24, 2015
 *      			Author: Minkwan Kim
 *
 * BoundaryConditionOther.cpp
 * 			-  
 *  
 */



#include "CFD/include/BoundaryConditions.hpp"
#include "CFD/include/Variables.hpp"
#include "Math/include/OP2A_Vector.hpp"

namespace OP2A{
namespace CFD{

void BCViscous::OtherTypeBC(Data::DataStorageVector<Data::DataStorage>& CellData1D_cl, Data::DataStorageVector<Data::DataStorage>& CellData1D_cr, int ND)
{
	CellData1D_cr	= CellData1D_cl;
}




}
}
