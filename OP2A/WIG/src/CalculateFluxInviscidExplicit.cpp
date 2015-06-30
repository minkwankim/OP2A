/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 29, 2015
 *      			Author: Minkwan Kim
 *
 * CalculateFluxInviscidImplicit.cpp
 * 			-  
 *  
 */


#include <omp.h>
#include "CFD/include/VariableSet.hpp"
#include "CFD/include/MiscFunctionsCFD.hpp"
#include "CFD/include/FluxInviscid.hpp"

#include "Math/include/MathMisc.hpp"

#include "../include/OP2A_Application.hpp"


void ApplicationOP2A::CalculateFluxInviscidExplicit()
{
	if (problem_setup.NUMERICAL_ORDER == 1)
	{
		CalculateFluxInviscidExplicit_SWFVM_1stOrder();
	}
	else
	{
		CalculateFluxInviscidExplicit_SWFVM_MUSCL();
	}
}



void ApplicationOP2A::CalculateFluxInviscidExplicit_SWFVM_MUSCL()
{

	int VAR	= grid.cells[1].data1D(indexW).numData;
	int indexE = species_set.NS + grid.ND;

#pragma omp parallel for num_threads(CFD_NT)
	for (int f = 1; f <= grid.NFM; f++)
	{
		Data::DataStorageVector<Data::DataStorage> data1D_L = cell_data1D_template;
		Data::DataStorageVector<Data::DataStorage> data1D_R = cell_data1D_template;

		CFD::Reconstruct::SecondOrderMUSCL(grid.faces[f].geo.cl[GRID::StencilLabel::CLL]->data1D(indexW),
											grid.faces[f].geo.cl[GRID::StencilLabel::CL]->data1D(indexW),
											grid.faces[f].geo.cr[GRID::StencilLabel::CR]->data1D(indexW),
											grid.faces[f].geo.cr[GRID::StencilLabel::CRR]->data1D(indexW),
											grid.faces[f].geo.cl[GRID::StencilLabel::CLL]->geo.x,
											grid.faces[f].geo.cl[GRID::StencilLabel::CL]->geo.x,
											grid.faces[f].geo.cr[GRID::StencilLabel::CR]->geo.x,
											grid.faces[f].geo.cr[GRID::StencilLabel::CRR]->geo.x,
											grid.faces[f].geo.x, true, problem_setup.LIMITER,
											data1D_L(indexW), data1D_R(indexW));

		CFD::VariableChange::W_to_V(CFD_variabletype, CFD_NT, data1D_L(indexW), species_set, grid.ND, data1D_L(indexV));
		CFD::VariableChange::V_to_Q(CFD_variabletype, CFD_NT, data1D_L(indexV), species_set, grid.ND, data1D_L(indexQ));

		CFD::VariableChange::W_to_V(CFD_variabletype, CFD_NT, data1D_L(indexW), species_set, grid.ND, data1D_R(indexV));
		CFD::VariableChange::V_to_Q(CFD_variabletype, CFD_NT, data1D_L(indexV), species_set, grid.ND, data1D_R(indexQ));


		double p_cl	= grid.faces[f].geo.cl[GRID::StencilLabel::CL]->data1D(indexW)(indexE);
		double p_cr	= grid.faces[f].geo.cr[GRID::StencilLabel::CR]->data1D(indexW)(indexE);
		double dp	= Math::fabs<double>(p_cl - p_cr) / Math::fmin<double>(p_cl, p_cr);

		CFD::FluxInviscid::SWFVS_Explicit(data1D_L, data1D_R, species_set, grid.ND,
										CFD_variabletype, indexQ, indexV, indexW,
										grid.faces[f].geo.n, f,
										dp, grid.faces[f].geo.dist_wall, grid.faces[f].geo.n_dot_wall, 5.0, 1.0e-5, 0.3,
										grid.faces[f].data1D(0));
	}
}


void ApplicationOP2A::CalculateFluxInviscidExplicit_SWFVM_1stOrder()
{

	int VAR	= grid.cells[1].data1D(indexW).numData;
	int indexE = species_set.NS + grid.ND;

#pragma omp parallel for
	for (int f = 1; f <= grid.NFM; f++)
	{
		Data::DataStorageVector<Data::DataStorage> data1D_L = cell_data1D_template;
		Data::DataStorageVector<Data::DataStorage> data1D_R = cell_data1D_template;

		CFD::Reconstruct::FirstOrder(grid.faces[f].geo.cl[GRID::StencilLabel::CL]->data1D(indexW),
									grid.faces[f].geo.cr[GRID::StencilLabel::CR]->data1D(indexW),
									data1D_L(indexW), data1D_R(indexW));

		CFD::VariableChange::W_to_V(CFD_variabletype, CFD_NT, data1D_L(indexW), species_set, grid.ND, data1D_L(indexV));
		CFD::VariableChange::V_to_Q(CFD_variabletype, CFD_NT, data1D_L(indexV), species_set, grid.ND, data1D_L(indexQ));

		CFD::VariableChange::W_to_V(CFD_variabletype, CFD_NT, data1D_R(indexW), species_set, grid.ND, data1D_R(indexV));
		CFD::VariableChange::V_to_Q(CFD_variabletype, CFD_NT, data1D_R(indexV), species_set, grid.ND, data1D_R(indexQ));

		double p_cl	= grid.faces[f].geo.cl[GRID::StencilLabel::CL]->data1D(indexW)(indexE);
		double p_cr	= grid.faces[f].geo.cr[GRID::StencilLabel::CR]->data1D(indexW)(indexE);
		double dp	= Math::fabs<double>(p_cl - p_cr) / Math::fmin<double>(p_cl, p_cr);


		CFD::FluxInviscid::SWFVS_Explicit(data1D_L, data1D_R, species_set, grid.ND,
										CFD_variabletype, indexQ, indexV, indexW,
										grid.faces[f].geo.n, f,
										dp, grid.faces[f].geo.dist_wall, grid.faces[f].geo.n_dot_wall, 5.0, 1.0e-5, 0.3,
										grid.faces[f].data1D(0));
	}
}
