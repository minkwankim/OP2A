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

#pragma omp parallel for
	for (int s = 0; 			s <= index_rhou-1; 		s++)	Qcr(s)	= Qcl(s);

#pragma omp parallel for
	for (int k = index_rhou;	k <= index_rhou+ND-1;	k++)	Qcr(k)	= Qcl(k);

#pragma omp parallel for
	for (int n = index_E;		n <= Qcl.data.size()-1;	n++)	Qcr(n)	= Qcl(n);
}







}
}
