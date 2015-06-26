/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 23, 2015
 *      			Author: Minkwan Kim
 *
 * ApplyBCInviscidNormal.cpp
 * 			-  
 *  
 */



#include "CFD/include/BoundaryConditions.hpp"
#include "CFD/include/Variables.hpp"

#include "../include/OP2A_Application.hpp"




void ApplicationOP2A::ApplyBCInviscidNormal()
{
	int indexAMB = 0;

#pragma omp parallel for num_threads(CFD_NT)
	for (int c = 1; c <= grid.NGM; c++)
	{
		CFD::CFDBCTypes BC_index = CFD::BCInviscid::BCTypeInCFD(grid.cells_ghost[c].geo.face_list[0]->geo.BC);

		switch (BC_index)
		{
		case CFD::CFDBCTypes::WallType:
			CFD::BCInviscid::wallTypeBC(grid.cells_ghost[c].geo.face_list[0]->geo.cl[0]->data1D(indexQ),
										grid.cells_ghost[c].geo.face_list[0]->geo.n,
										grid.cells_ghost[c].geo.face_list[0]->geo.cr[0]->data1D(indexQ), grid.ND);
			break;

		case CFD::CFDBCTypes::InletType:
			CFD::BCInviscid::inletTypeBC(grid.cells_ghost[c].geo.face_list[0]->geo.cl[0]->data1D(indexQ),
										grid.cells_ghost[c].geo.face_list[0]->geo.cr[0]->data1D(indexQ),
										grid.ND,
										IC_Q(1));
			break;

		case CFD::CFDBCTypes::FreestreamType:
			CFD::BCInviscid::inletTypeBC(grid.cells_ghost[c].geo.face_list[0]->geo.cl[0]->data1D(indexQ),
										grid.cells_ghost[c].geo.face_list[0]->geo.cr[0]->data1D(indexQ),
										grid.ND,
										IC_Q(indexAMB));
			break;

		case CFD::CFDBCTypes::ExitType:
			CFD::BCInviscid::exitTypeBC(grid.cells_ghost[c].geo.face_list[0]->geo.cl[0]->data1D(indexQ),
										grid.cells_ghost[c].geo.face_list[0]->geo.cr[0]->data1D(indexQ),
										grid.ND);
			break;
		}
	}


#pragma omp parallel for num_threads(CFD_NT)
	for (int c = 1; c <= grid.NGM; c++)
	{
		CFD::CellData1D::CompleteDataWIGCase1(grid.cells_ghost[c].data1D, species_set, grid.ND, CFD_variabletype, CFD_NT, indexQ, indexV, indexW, indexMIX, indexXs, indexYs);
	}
}

