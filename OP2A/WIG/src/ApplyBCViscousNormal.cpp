/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 24, 2015
 *      			Author: Minkwan Kim
 *
 * ApplyBCViscousNorman.cpp
 * 			-  
 *  
 */



#include "CFD/include/BoundaryConditions.hpp"
#include "CFD/include/Variables.hpp"
#include "CFD/include/VariableSet.hpp"

#include "../include/OP2A_Application.hpp"



void ApplicationOP2A::ApplyBCViscousNormal()
{
	int indexAMB = 0;
	int indexkappa_tra = grid.cells[1].data1D(indexMIX).dataMap.find(CFD::CFD_VariableSet::kappa_mix(CHEM::EnergyMode::E_TRA));
	int NE	= CFD::VariableChange::NumEnergyModes(CFD_variabletype);

#pragma omp parallel for num_threads(CFD_NT)
	for (int c = 1; c <= grid.NGM; c++)
	{
		unsigned int wallID;

		CFD::CFDBCTypes BC_indexViscous = CFD::BCViscous::BCTypeInCFD(grid.cells_ghost[c].geo.face_list[0]->geo.BC);

		if  (BC_indexViscous == CFD::CFDBCTypes::WallType)
		{
			wallID = 0;
			double kappa = 0.0;

#pragma omp parallel for reduction(+:kappa)
			for (int em = indexkappa_tra; em <= indexkappa_tra+NE-1; em++)
			{
				kappa	+= grid.cells_ghost[c].geo.face_list[0]->geo.cl[0]->data1D(indexMIX)(em);
			}


			CFD::BCViscous::wallTypeBC(grid.cells_ghost[c].geo.face_list[0]->geo.cl[0]->data1D,
										grid.cells_ghost[c].geo.face_list[0]->geo.cr[0]->data1D,
										species_set,
										BCValuesWall(wallID),
										grid.cells_ghost[c].geo.face_list[0]->geo.cl[0]->geo.x,
										grid.cells_ghost[c].geo.face_list[0]->geo.x,
										grid.ND,
										CFD_variabletype,
										CFD_NT,
										problem_setup.adiabatic,
										problem_setup.catalytic,
										true,
										problem_setup.radiative,
										kappa,
										problem_setup.use_emissivity,
										problem_setup.emissivity_Tref[wallID],
										problem_setup.emissivity_below[wallID],
										problem_setup.emissivity_above[wallID]);
		}
		else if (BC_indexViscous == CFD::CFDBCTypes::OtherType)
		{

			CFD::BCViscous::OtherTypeBC(grid.cells_ghost[c].geo.face_list[0]->geo.cl[0]->data1D,
										grid.cells_ghost[c].geo.face_list[0]->geo.cr[0]->data1D,
										grid.ND);
		}
	}

}
