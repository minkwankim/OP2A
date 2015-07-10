/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 23, 2015
 *      			Author: Minkwan Kim
 *
 * 			-  
 *  
 */

#include "Common/include/MultiDimension.hpp"
#include "CFD/include/BoundaryConditions.hpp"

namespace OP2A{
namespace CFD{

void  BCInviscid::wallTypeBC(Data::DataStorage& Qcl, vector< vector<double> > & face_normal_vectors, Data::DataStorage& Qcr, int ND)
{
	vector<double>	rhoV(ND, 0.0);
	vector<double>	rhoVp(ND, 0.0);

	int index_rhou	= Qcl.dataMap.find(CFD_VariableSet::var_stringQ(0, FluxCategory::Momentum));
	int index_E		= Qcl.dataMap.find(CFD_VariableSet::var_stringQ(CHEM::EnergyMode::E_TRA, FluxCategory::Energy));


	for (int k1 = 0; k1 <= ND-1; k1++)
	{
		for (int k2 = 0; k2 <= ND-1; k2++)
		{
			rhoVp[k1]	+= Qcl(index_rhou+k2) * face_normal_vectors[k1][k2];
		}
	}
	rhoVp[0]	=	-rhoVp[0];


	for (int k2 = 0; k2 <= ND-1; k2++)
	{
		for (int k1 = 0; k1 <= ND-1; k1++)
		{
			rhoV[k2]	+= rhoVp[k1] * face_normal_vectors[k1][k2];
		}
	}

//#pragma omp parallel for
#pragma ivdep
	for (int s = 0; 			s <= index_rhou-1; 		s++)	Qcr(s)	= Qcl(s);

//#pragma omp parallel for
#pragma ivdep
	for (int k = 0;	k <= ND-1;	k++)	Qcr(index_rhou+k)	= rhoV[k];

//#pragma omp parallel for
#pragma ivdep
	for (int n = index_E;		n <= Qcl.data.size()-1;	n++)	Qcr(n)	= Qcl(n);
}



void  BCInviscidImplicit::wallTypeBC(Data::DataStorage2D& J_plus, Data::DataStorage2D& J_minus, vector< vector<double> >& face_normal_vector, int NS, int ND, int NE)
{
	int VAR = NS + ND + NE;

	vector< vector<double> >	E 		= Common::vector_2D(VAR, VAR, 0.0);
	vector< vector<double> >	R 		= Common::vector_2D(ND, ND, 0.0);
	vector< vector<double> >	R_inv 	= Common::vector_2D(ND, ND, 0.0);

#pragma ivdep
	for (int i = 0; i <= VAR-1; i++)
	{
		E[i][i]	= 1.0;
	}


#pragma ivdep
	for (int j = 0; j <= ND-1; j++)
	{
		for (int k = 0; k <= ND-1; k++)
		{
			R[j][k]		= face_normal_vector[j][k];
			R_inv[j][k]	= face_normal_vector[k][j];
		}
	}

#pragma ivdep
	for (int k = 0; k <= ND-1; k++)
	{
		R[0][k]	= -R[0][k];
	}

#pragma ivdep
	for (int k = 0; k <= ND-1; k++)
	{
		for (int l = 0; l <= ND-1; l++)
		{
			double aux = 0.0;
			for (int m = 0; m <= ND-1; m++)
			{
				aux	+= R_inv[k][m] * R[m][l];
			}

			E[NS + k][NS + l]	= aux;
		}
	}

#pragma ivdep
	for (int k = 0; k <= VAR-1; k++)
	{
		for (int l = 0; l <= VAR-1; l++)
		{
			double aux = 0.0;
			for (int m = 0; m <= VAR-1; m++)
			{
				aux	+= J_minus(k,m) * E[m][l];
			}

			J_plus(k,l)		+= aux;
			J_minus(k, l)	= 0.0;
		}
	}
}




}
}
