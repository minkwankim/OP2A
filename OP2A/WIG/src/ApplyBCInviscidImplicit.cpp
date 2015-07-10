/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 9, 2015
 *      			Author: Minkwan Kim
 *
 * ApplyBCInviscidImplicit.cpp
 * 			-  
 *  
 */




#include "CFD/include/BoundaryConditions.hpp"
#include "CFD/include/Variables.hpp"

#include "../include/OP2A_Application.hpp"




void ApplicationOP2A::ApplyBCInviscidImplicit()
{
	int NS  = species_set.NS;
	int ND	= grid.ND;
	int NE  = grid.cells[1].data1D(0).numData - NS - ND;
	int index_J_plus  = 0;
	int index_J_minus = 1;


#pragma omp parallel for num_threads(CFD_NT)
	for (int c = 1; c <= grid.NGM; c++)
	{
		CFD::CFDBCTypes BC_index = CFD::BCInviscid::BCTypeInCFD(grid.cells_ghost[c].geo.face_list[0]->geo.BC);

		switch (BC_index)
		{
		case CFD::CFDBCTypes::WallType:
			CFD::BCInviscidImplicit::wallTypeBC(grid.cells_ghost[c].geo.face_list[0]->data2D(index_J_plus),
												grid.cells_ghost[c].geo.face_list[0]->data2D(index_J_minus),
												grid.cells_ghost[c].geo.face_list[0]->geo.n,
												NS, ND, NE);
			break;

		case CFD::CFDBCTypes::InletType:
			CFD::BCInviscidImplicit::inletTypeBC(grid.cells_ghost[c].geo.face_list[0]->data2D(index_J_plus),
												grid.cells_ghost[c].geo.face_list[0]->data2D(index_J_minus),
												grid.cells_ghost[c].geo.face_list[0]->geo.n,
												NS, ND, NE);
			break;

		case CFD::CFDBCTypes::FreestreamType:
			CFD::BCInviscidImplicit::inletTypeBC(grid.cells_ghost[c].geo.face_list[0]->data2D(index_J_plus),
												grid.cells_ghost[c].geo.face_list[0]->data2D(index_J_minus),
												grid.cells_ghost[c].geo.face_list[0]->geo.n,
												NS, ND, NE);
			break;

		case CFD::CFDBCTypes::ExitType:
			CFD::BCInviscidImplicit::exitTypeBC(grid.cells_ghost[c].geo.face_list[0]->data2D(index_J_plus),
												grid.cells_ghost[c].geo.face_list[0]->data2D(index_J_minus),
												grid.cells_ghost[c].geo.face_list[0]->geo.n,
												NS, ND, NE);
			break;
		}
	}
}
