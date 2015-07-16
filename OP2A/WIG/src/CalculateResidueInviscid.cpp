/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 29, 2015
 *      			Author: Minkwan Kim
 *
 * CalculateResidueInviscid.cpp
 * 			-  
 *  
 */



#include <omp.h>
#include <limits>
#include "Common/include/Exception_NaNValue.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"

#include "Math/include/MathMisc.hpp"

#include "../include/OP2A_Application.hpp"



void ApplicationOP2A::CalculateResidueInviscid()
{
	if (omp_get_num_threads() > 4)
	{
		CalculateResidueInviscid_ver1();
	}
	else
	{
		CalculateResidueInviscid_ver1();
	}
}



void ApplicationOP2A::CalculateResidueInviscid_ver1()
{


	// Initialize
#pragma omp parallel num_threads(CFD_NT)
	for (int c = 1; c <= grid.NCM; c++)
	{
#pragma ivdep
		for (int i = 0; i <= grid.cells[c].data1D(indexResidue).numData-1; i++)
		{
			grid.cells[c].data1D(indexResidue)(i)	= 0.0;
		}
	}



	// Calculate Residue
	if (problem_setup.is_axisymmetric == true)
	{
#pragma omp parallel for num_threads(CFD_NT)
		for (int c = 1; c <= grid.NCM; c++)
		{
			for (int index_i = 0; index_i <= grid.cells[c].geo.NF-1; index_i++)
			{
				double S = 0.0;
				S = grid.cells[c].geo.face_list[index_i]->geo.S * Math::fabs<double>(grid.cells[c].geo.face_list[index_i]->geo.x[1]);

				for (int k = 0; k <= grid.cells[c].geo.face_list[index_i]->data1D(0).numData-1; k++)
				{
					double flux	= grid.cells[c].geo.face_list[index_i]->data1D(0)(k) * S;

					if (grid.cells[c].geo.face_list[index_i]->geo.cr[0]->geo.ID == grid.cells[c].geo.ID)
					{
						grid.cells[c].data1D(indexResidue)(k)	-= flux;
					}
					else if (grid.cells[c].geo.face_list[index_i]->geo.cl[0]->geo.ID == grid.cells[c].geo.ID)
					{
						grid.cells[c].data1D(indexResidue)(k)	+= flux;
					}
				}
			}
		}


		int indexAxis = species_set.NS +1;
		int indexE = species_set.NS +grid.ND;

#pragma omp parallel for num_threads(CFD_NT)
		for (int c = 1;  c <= grid.NCM; c++)
		{
			double	S	= grid.cells[c].geo.S;
			grid.cells[c].data1D(indexResidue)(indexAxis)	-= S * grid.cells[c].data1D(indexW)(indexE) * Math::fsgn<double>(grid.cells[c].geo.x[1]);
		}
	}
	else
	{
#pragma omp parallel for num_threads(CFD_NT)
		for (int c = 1; c <= grid.NCM; c++)
		{
			for (int index_i = 0; index_i <= grid.cells[c].geo.NF-1; index_i++)
			{
				double S = 0.0;
				S = grid.cells[c].geo.face_list[index_i]->geo.S;

				for (int k = 0; k <= grid.cells[c].geo.face_list[index_i]->data1D(0).numData-1; k++)
				{
					double flux	= grid.cells[c].geo.face_list[index_i]->data1D(0)(k) * S;

					if (grid.cells[c].geo.face_list[index_i]->geo.cr[0]->geo.ID == grid.cells[c].geo.ID)
					{
						grid.cells[c].data1D(indexResidue)(k)	-= flux;
					}
					else if (grid.cells[c].geo.face_list[index_i]->geo.cl[0]->geo.ID == grid.cells[c].geo.ID)
					{
						grid.cells[c].data1D(indexResidue)(k)	+= flux;
					}
				}
			}
		}
	}


	// Check Error
#pragma omp parallel for num_threads(CFD_NT)
	for (int c = 1; c <= grid.NCM; c++)
	{

#pragma ivdep
		for (int k = 0; k <= grid.cells[c].data1D(indexResidue).numData-1; k++)
		{

			/*
			if (Math::fabs<double>(grid.cells[c].data1D(indexResidue)(k)) <= 1.0e-18)
			{
				grid.cells[c].data1D(indexResidue)(k) = 0.0;
			}
			*/

			if (grid.cells[c].data1D(indexResidue)(k) != grid.cells[c].data1D(indexResidue)(k))
			{
				std::ostringstream oss;
				oss << "NaN value in the calculation of inviscid Residue: [Cell ID]: " << c << "  [VAR]: " << k;
				throw Common::ExceptionNaNValue (FromHere(), oss.str());
			}

			if (grid.cells[c].data1D(indexResidue)(k) == numeric_limits<double>::infinity())
			{
				std::ostringstream oss;
				oss << "Infinite Value in the calculation of inviscid Residue: [Cell ID]: " << c << "  [VAR]: " << k;
				throw Common::ExceptionInfiniteValue (FromHere(), oss.str());
			}

		}
	}

}







void ApplicationOP2A::CalculateResidueInviscid_ver2()
{


	// Initialize
#pragma omp parallel
	for (int c = 1; c <= grid.NCM; c++)
	{
#pragma ivdep
		for (int i = 0; i <= grid.cells[c].data1D(indexResidue).numData-1; i++)
		{
			grid.cells[c].data1D(indexResidue)(i)	= 0.0;
		}
	}

#pragma omp parallel
	for (int c = 1; c <= grid.NGM; c++)
	{
#pragma ivdep
		for (int i = 0; i <= grid.cells_ghost[c].data1D(indexResidue).numData-1; i++)
		{
			grid.cells_ghost[c].data1D(indexResidue)(i)	= 0.0;
		}
	}


	// Calculate Residue
	if (problem_setup.is_axisymmetric == true)
	{
		for (int f = 1; f <= grid.NFM; f++)
		{
			double S = 0.0;
			S = grid.faces[f].geo.S * Math::fabs<double>(grid.faces[f].geo.x[1]);

#pragma ivdep
			for (int k = 0; k <= grid.faces[f].data1D(0).numData-1; k++)
			{
				double flux	= grid.faces[f].data1D(0)(k) * S;

				grid.faces[f].geo.cl[0]->data1D(indexResidue)(k)	+= flux;
				grid.faces[f].geo.cr[0]->data1D(indexResidue)(k)	-= flux;
			}
		}

		int indexAxis = species_set.NS +1;
		int indexE = species_set.NS +grid.ND;

#pragma omp parallel for
		for (int c = 1;  c <= grid.NCM; c++)
		{
			double	S	= grid.cells[c].geo.S;
			grid.cells[c].data1D(indexResidue)(indexAxis)	-= S * grid.cells[c].data1D(indexW)(indexE) * Math::fsgn<double>(grid.cells[c].geo.x[1]);
		}
	}
	else
	{
		for (int f = 1; f <= grid.NFM; f++)
		{
			double S = 0.0;
			S = grid.faces[f].geo.S;;

#pragma ivdep
			for (int k = 0; k <= grid.faces[f].data1D(0).numData-1; k++)
			{
				double flux	= grid.faces[f].data1D(0)(k) * S;

				grid.faces[f].geo.cl[0]->data1D(indexResidue)(k)	+= flux;
				grid.faces[f].geo.cr[0]->data1D(indexResidue)(k)	-= flux;
			}
		}
	}


	// Check Error
#pragma omp parallel for
	for (int c = 1; c <= grid.NCM; c++)
	{
#pragma ivdep
		for (int k = 0; k <= grid.cells[c].data1D(indexResidue).numData-1; k++)
		{

			if (grid.cells[c].data1D(indexResidue)(k) != grid.cells[c].data1D(indexResidue)(k))
			{
				std::ostringstream oss;
				oss << "NaN value in the calculation of inviscid Residue: [Cell ID]: " << c << "  [VAR]: " << k;
				throw Common::ExceptionNaNValue (FromHere(), oss.str());
			}

			if (grid.cells[c].data1D(indexResidue)(k) == numeric_limits<double>::infinity())
			{
				std::ostringstream oss;
				oss << "Infinite Value in the calculation of inviscid Residue: [Cell ID]: " << c << "  [VAR]: " << k;
				throw Common::ExceptionInfiniteValue (FromHere(), oss.str());
			}

		}
	}

}
