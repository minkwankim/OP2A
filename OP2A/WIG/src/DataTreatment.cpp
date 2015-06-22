/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 22, 2015
 *      			Author: Minkwan Kim
 *
 * CellDataTreatment.cpp
 * 			-  
 *  
 */



#include "Common/include/CodeLocation.hpp"
#include "Common/include/Exception_NPExceed.hpp"
#include "Common/include/StringOps.hpp"
#include "CFD/include/VariableSamplesWIG.hpp"

#include "CFD/include/Variables.hpp"
#include "Math/include/MathMisc.hpp"

#include "../include/OP2A_Application.hpp"



void ApplicationOP2A::Cell1DDataTreatement(int typeCase, bool is_initialize)
{
	int indexQ		= grid.cells[1].data1D.dataMap.find(NAME_Q);
	int indexV		= grid.cells[1].data1D.dataMap.find(NAME_V);
	int indexW		= grid.cells[1].data1D.dataMap.find(NAME_W);
	int indexMIX	= grid.cells[1].data1D.dataMap.find(NAME_MIX);
	int indexXs		= grid.cells[1].data1D.dataMap.find(NAME_XS);
	int indexYs		= grid.cells[1].data1D.dataMap.find(NAME_YS);
	int indexR		= grid.cells[1].data1D.dataMap.find(NAME_R);
	int indexdQ		= grid.cells[1].data1D.dataMap.find(NAME_dQ);
	int indexQnew	= grid.cells[1].data1D.dataMap.find(NAME_Qnew);
	int indexS		= grid.cells[1].data1D.dataMap.find(NAME_S);

	int indexDivV	= -1;
	int indexdpdQ	= -1;
	int indexDs		= -1;
	int indexMus	= -1;


	if (problem_setup.is_axisymmetric == true) 		indexDivV	= grid.cells[1].data1D.dataMap.find(NAME_DIVV);
	if (problem_setup.TIME_INTEGRATION_METHOD != 0)	indexdpdQ	= grid.cells[1].data1D.dataMap.find(NAME_dpdQ);
	if (problem_setup.is_viscous == true)
	{
		indexDs		= grid.cells[1].data1D.dataMap.find(NAME_DS);
		indexMus	= grid.cells[1].data1D.dataMap.find(NAME_MUS);
	}


#pragma omp parallel for num_threads(CFD_NT)
	for (int c = 1; c <= grid.NCM; c++)
	{
		CFD::CellData1D::CompleteDataWIG(grid.cells[c].data1D, species_set, grid.ND, CFD_variabletype, CFD_NT, indexQ, indexV, indexW, indexMIX, indexXs, indexYs, indexR, indexdQ, indexQnew, indexS, indexDivV, indexdpdQ, indexDs,indexMus, typeCase, is_initialize);
	}
}


void ApplicationOP2A::Cell2DDataTreatement(int typeCase, bool is_initialize)
{
	int indexdTdQ	= -1;
	int indexdSdQ	= -1;
	int indexKappas	= -1;


	if (problem_setup.is_viscous == true)	indexKappas	= grid.cells[1].data2D.dataMap.find(NAME_KAPPAS);
	if (problem_setup.TIME_INTEGRATION_METHOD != 0)
	{
		indexdTdQ	= grid.cells[1].data2D.dataMap.find(NAME_dTdQ);
		indexdSdQ	= grid.cells[1].data2D.dataMap.find(NAME_dSdQ);
	}


#pragma omp parallel for num_threads(CFD_NT)
	for (int c = 1; c <= grid.NCM; c++)
	{
		CFD::CellData2D::InitializeData(grid.cells[c].data2D, species_set, grid.ND, CFD_variabletype, CFD_NT, indexdTdQ, indexdSdQ, indexKappas);
	}
}


void ApplicationOP2A::Face1DDataTreatement(int typeCase, bool is_initialize)
{
#pragma omp parallel for num_threads(CFD_NT)
	for (int f = 1; f <= grid.NFM; f++)
	{
		CFD::FaceData1D::InitializeData(grid.faces[f].data1D, species_set, grid.ND, CFD_variabletype, CFD_NT, problem_setup.is_viscous);
	}
}


void ApplicationOP2A::Face2DDataTreatement(int typeCase, bool is_initialize)
{
#pragma omp parallel for num_threads(CFD_NT)
	for (int f = 1; f <= grid.NFM; f++)
	{
		CFD::FaceData2D::InitializeData(grid.faces[f].data2D, species_set, grid.ND, CFD_variabletype, CFD_NT, problem_setup.is_viscous, problem_setup.TIME_INTEGRATION_METHOD);
	}
}



void ApplicationOP2A::DataTreatement(int typeCase, bool is_initialize, int type)
{
	switch (type)
	{
	case OP2A_DataCatogories::OP2A_CELL1D:
		Cell1DDataTreatement(typeCase, is_initialize);
		break;

	case OP2A_DataCatogories::OP2A_CELL2D:
		Cell2DDataTreatement(typeCase, is_initialize);
		break;

	case OP2A_DataCatogories::OP2A_FACE1D:
		Face1DDataTreatement(typeCase, is_initialize);
		break;

	case OP2A_DataCatogories::OP2A_FACE2D:
		Face2DDataTreatement(typeCase, is_initialize);
		break;
	}

}

