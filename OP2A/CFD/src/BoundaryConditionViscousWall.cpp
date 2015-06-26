/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 23, 2015
 *      			Author: Minkwan Kim
 *
 * BoundaryConditionViscousWall.cpp
 * 			-  
 *  
 */



#include "CFD/include/BoundaryConditions.hpp"
#include "CFD/include/Variables.hpp"
#include "Math/include/OP2A_Vector.hpp"

namespace OP2A{
namespace CFD{


void  BCViscous::wallTypeBC(Data::DataStorageVector<Data::DataStorage>& CellData1D_cl,
							Data::DataStorageVector<Data::DataStorage>& CellData1D_cr,
							CHEM::SpeciesSet& species_set,
							Data::DataStorage& BCValuesWall,
							vector<double>& xcl,
							vector<double>& xf,
							int  ND,
							int  CFD_variabletype,
							int	 CFD_NT,
							bool adiabaticWall,
							bool catalyticWall,
							bool nonSlipWall,
							bool radiativeWall,
							double kappa_cl,
							bool use_emissivity,
							double Tref,
							vector<double> emissivity_below,
							vector<double> emissivity_above)
{

	int iterMax = 5000;
	double eps0 = 0.001;

	// Getting initial information
	int indexQ 		= CellData1D_cl.dataMap.find(NAME_Q);
	int indexV 		= CellData1D_cl.dataMap.find(NAME_V);
	int indexW 		= CellData1D_cl.dataMap.find(NAME_W);
	int indexMIX 	= CellData1D_cl.dataMap.find(NAME_MIX);
	int indexXs 	= CellData1D_cl.dataMap.find(NAME_XS);
	int indexYs 	= CellData1D_cl.dataMap.find(NAME_YS);

	int index_u		= CellData1D_cl(indexV).dataMap.find(CFD_VariableSet::var_stringV(0, FluxCategory::Momentum));
	int index_T		= CellData1D_cl(indexV).dataMap.find(CFD_VariableSet::var_stringV(CHEM::EnergyMode::E_TRA, FluxCategory::Energy));


	// Case 1: Adiabatic wall
	if (adiabaticWall == true)
	{
#pragma omp parallel for
		for (int s = 0;	s <= index_u-1; s++)	CellData1D_cr(indexV)(s)	= CellData1D_cl(indexV)(s);

#pragma omp parallel for
		for (int k = index_u; k <= index_u+ND-1; k++) CellData1D_cr(indexV)(k)	= -CellData1D_cl(indexV)(k);

#pragma omp parallel for
		for (int n = index_T; n <= CellData1D_cl(indexV).numData-1; n++) CellData1D_cr(indexV)(n)	= CellData1D_cl(indexV)(n);

		CellData1D::CompleteDataWIGCase2(CellData1D_cr, species_set, ND, CFD_variabletype, CFD_NT, indexQ, indexV, indexW, indexMIX, indexXs, indexYs);
	}
	// Non-Adiabatic wall
	else
	{
		Data::DataStorage Ys_wall("Ys at Wall", CellData1D_cl(indexYs).numData);					// Ys at wall
		Data::DataStorage V_wall("Primitive value at wall face", CellData1D_cl(indexV).numData);	// V at  wall
		double p_wall;																				// pressure at wall

		// A. Mass Fraction
		if (catalyticWall == true)
		{
#pragma omp parallel for
			for (int s = 0; s <= Ys_wall.numData-1; s++)	Ys_wall(s)	= BCValuesWall(s);
		}
		else
		{
#pragma omp parallel for
			for (int s = 0; s <= Ys_wall.numData-1; s++)	Ys_wall(s)	= CellData1D_cl(indexYs)(s);
		}


		// B. NonSlip Wall
		if (nonSlipWall)
		{
			// 1. Velocities
#pragma omp parallel for
			for (int k = index_u;	k <= index_u +ND-1;	k++)	V_wall(k)	= 0.0;

			// 2. Temperatures
			if (radiativeWall == true)
			{
				Math::VECTOR x_f_cl(xf, xcl);
				double dn = x_f_cl.length();

				V_wall(index_T)	= CFD_CalculateWallTemperature(CellData1D_cl(indexV)(index_T), kappa_cl, dn, use_emissivity, Tref, emissivity_below, emissivity_above, iterMax, eps0);

				#pragma omp parallel for
				for (int n = index_T+1; n <= CellData1D_cl(indexV).numData-1; n++) V_wall(n)	= V_wall(index_T);
			}
			else
			{
#pragma omp parallel for
				for (int n = index_T; n <= CellData1D_cl(indexV).numData-1; n++) V_wall(n)	= BCValuesWall(n);
			}
		}
		else
		{
			/*
			 * @todo	Need to add slip wall BC
			 */
		}


		// 1. Find Values at cr
		//	- pressure
		double p_cr	= CellData1D_cl(indexW)(index_T);

		// 	- Mass Fraction
#pragma omp parallel for
		for (int s = 0;	s <= index_u-1; s++)	CellData1D_cr(indexYs)(s)	= CellData1D_cl(indexYs)(s) - Ys_wall(s);

		//	- Velocity
#pragma omp parallel for
		for (int k = index_u; k <= index_u+ND-1; k++) CellData1D_cr(indexV)(k)	= 2.0*V_wall(k) - CellData1D_cl(indexV)(k);

		//	- Temperature
#pragma omp parallel for
		for (int n = index_T; n <= CellData1D_cl(indexV).numData-1; n++)
		{
			double Ttemp	= 2.0*V_wall(n) - CellData1D_cl(indexV)(n);

			if (Ttemp > 0.0)	CellData1D_cr(indexV)(n)	= Ttemp;
			else				CellData1D_cr(indexV)(n)	= V_wall(n);
		}

		double rho_cr = 0.0;
		double YsRsT = 0.0;
		double Tcr, Tecr;

		switch (CFD_variabletype)
		{
		case 1:
			Tcr	= CellData1D_cr(indexV)(index_T);

#pragma omp parallel for reduction(+:YsRsT)
			for (int s = 0;	s <= index_u-1; s++)
			{
				YsRsT	+= CellData1D_cl(indexYs)(s) * species_set.species[s].R * Tcr;
			}

			rho_cr	= p_cr / YsRsT;
			break;

		case 2:
			Tcr	 = CellData1D_cr(indexV)(index_T);
			Tecr = CellData1D_cr(indexV)(index_T+1);
#pragma omp parallel for reduction(+:YsRsT)
			for (int s = 0;	s <= index_u-1; s++)
			{
				if (species_set.species[s].type != CHEM::SpeciesType::Electron)
				{
					YsRsT	+= CellData1D_cl(indexYs)(s) * species_set.species[s].R * Tcr;
				}
				else
				{
					YsRsT	+= CellData1D_cl(indexYs)(s) * species_set.species[s].R * Tecr;
				}
			}

			rho_cr	= p_cr / YsRsT;
			break;

		case 3:
			Tcr	= CellData1D_cr(indexV)(index_T);

#pragma omp parallel for reduction(+:YsRsT)
			for (int s = 0;	s <= index_u-1; s++)
			{
				YsRsT	+= CellData1D_cl(indexYs)(s) * species_set.species[s].R * Tcr;
			}

			rho_cr	= p_cr / YsRsT;
			break;

		case 4:
			Tcr	 = CellData1D_cr(indexV)(index_T);
			Tecr = CellData1D_cr(indexV)(index_T+2);
#pragma omp parallel for reduction(+:YsRsT)
			for (int s = 0;	s <= index_u-1; s++)
			{
				if (species_set.species[s].type != CHEM::SpeciesType::Electron)
				{
					YsRsT	+= CellData1D_cl(indexYs)(s) * species_set.species[s].R * Tcr;
				}
				else
				{
					YsRsT	+= CellData1D_cl(indexYs)(s) * species_set.species[s].R * Tecr;
				}
			}

			rho_cr	= p_cr / YsRsT;
			break;

		case 5:
			Tcr	= CellData1D_cr(indexV)(index_T);

#pragma omp parallel for reduction(+:YsRsT)
			for (int s = 0;	s <= index_u-1; s++)
			{
				YsRsT	+= CellData1D_cl(indexYs)(s) * species_set.species[s].R * Tcr;
			}

			rho_cr	= p_cr / YsRsT;
			break;

		case 6:
			Tcr	 = CellData1D_cr(indexV)(index_T);
			Tecr = CellData1D_cr(indexV)(index_T+2);
#pragma omp parallel for reduction(+:YsRsT)
			for (int s = 0;	s <= index_u-1; s++)
			{
				if (species_set.species[s].type != CHEM::SpeciesType::Electron)
				{
					YsRsT	+= CellData1D_cl(indexYs)(s) * species_set.species[s].R * Tcr;
				}
				else
				{
					YsRsT	+= CellData1D_cl(indexYs)(s) * species_set.species[s].R * Tecr;
				}
			}

			rho_cr	= p_cr / YsRsT;
			break;

		case 7:
			Tcr	= CellData1D_cr(indexV)(index_T);

#pragma omp parallel for reduction(+:YsRsT)
			for (int s = 0;	s <= index_u-1; s++)
			{
				YsRsT	+= CellData1D_cl(indexYs)(s) * species_set.species[s].R * Tcr;
			}

			rho_cr	= p_cr / YsRsT;
			break;

		case 8:
			Tcr	 = CellData1D_cr(indexV)(index_T);
			Tecr = CellData1D_cr(indexV)(index_T+3);
#pragma omp parallel for reduction(+:YsRsT)
			for (int s = 0;	s <= index_u-1; s++)
			{
				if (species_set.species[s].type != CHEM::SpeciesType::Electron)
				{
					YsRsT	+= CellData1D_cl(indexYs)(s) * species_set.species[s].R * Tcr;
				}
				else
				{
					YsRsT	+= CellData1D_cl(indexYs)(s) * species_set.species[s].R * Tecr;
				}
			}

			rho_cr	= p_cr / YsRsT;
			break;
		}


#pragma omp parallel for
		for (int s = 0;	s <= index_u-1; s++)	CellData1D_cr(indexV)(s)	= rho_cr * CellData1D_cr(indexYs)(s);

		CellData1D::CompleteDataWIGCase2(CellData1D_cr, species_set, ND, CFD_variabletype, CFD_NT, indexQ, indexV, indexW, indexMIX, indexXs, indexYs);
	}
}









}
}
