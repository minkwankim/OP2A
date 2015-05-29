/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * Solution_CFD.cpp
 * 			-  
 *  
 */



#include "../include/OP2A_DATA_CFD_solution.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"

SOL_CFD::SOL_CFD()
{

}


SOL_CFD::~SOL_CFD()
{

}



void SOL_CFD::allocate_size(int NNM, int NFM, int NCM, int NGM, bool viscous, int time_implicit, bool axisymmetric)
{
	Qc		= vector_2D(NCM, setup.VAR, 0.0);
	Vc		= vector_2D(NCM, setup.VAR, 0.0);
	Wc		= vector_2D(NCM, setup.VAR, 0.0);
	Rn		= vector_2D(NCM, setup.VAR, 0.0);
	S_source= vector_2D(NCM, setup.VAR, 0.0);
	dQ		= vector_2D(NCM, setup.VAR, 0.0);
	Q_new	= vector_2D(NCM, setup.VAR, 0.0);

	Qgc		= vector_2D(NGM, setup.VAR, 0.0);
	Vgc		= vector_2D(NGM, setup.VAR, 0.0);;
	Wgc		= vector_2D(NGM, setup.VAR, 0.0);;

	mixture_data_c.resize(NCM);
	mixture_data_gc.resize(NGM);

	Flux_f_inviscid = vector_2D(NFM, setup.VAR, 0.0);
	if (viscous == true)
	{
		Flux_f_viscous = vector_2D(NFM, setup.VAR, 0.0);

		transport_prop_c.resize(NCM);
		transport_prop_gc.resize(NGM);

		if (axisymmetric == true)	div_Vc.resize(NCM);
	}

	Preconditioner.resize(NCM);
	for (int c = 0; c <= NCM-1; c++)	Preconditioner[c]	=  vector_2D(setup.VAR, setup.VAR, 0.0);

	if (time_implicit != 0)
	{
		Jacobian_inviscid_plus.resize(NFM);
		Jacobian_inviscid_minus.resize(NFM);
		Jacobian_source.resize(NFM);

		dT_dQ.resize(NCM);
		dp_dQ.resize(NCM);
		for (int c = 0; c <= NCM-1; c++)
		{
			dT_dQ[c] = vector_2D(setup.NE, setup.VAR, 0.0);
			dp_dQ[c].resize(setup.VAR);
		}


		for (int f = 0; f <= NFM-1; f++)
		{
			Jacobian_inviscid_plus[f]	=  vector_2D(setup.VAR, setup.VAR, 0.0);
			Jacobian_inviscid_minus[f]	=  vector_2D(setup.VAR, setup.VAR, 0.0);
			Jacobian_source[f]			=  vector_2D(setup.VAR, setup.VAR, 0.0);
		}

		if (viscous == true)
		{
			Jacobian_viscous_plus.resize(NFM);
			Jacobian_viscous_minus.resize(NFM);

			for (int f = 0; f <= NCM-1; f++)
			{
				Jacobian_viscous_plus[f]	=  vector_2D(setup.VAR, setup.VAR, 0.0);
				Jacobian_viscous_minus[f]	=  vector_2D(setup.VAR, setup.VAR, 0.0);
			}
		}
	}
}
