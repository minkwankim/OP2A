/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 26, 2015
 *      			Author: Minkwan Kim
 *
 * Derivatives.cpp
 * 			-  
 *  
 */



#include <omp.h>

#include "CFD/include/Derivatives.hpp"


namespace OP2A{
namespace CFD{



void Derivatives::dTdQ(Data::DataStorageVector<Data::DataStorage>& data1D, CHEM::SpeciesSet& species_set, int ND,
				unsigned int type, unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX,
				Data::DataStorageVector<Data::DataStorage>& dT)
{
	int indexT  = 0;
	int indexTr = 0;
	int indexTv = 0;
	int indexTe = 0;

	switch (type)
	{
	case 1:
		DerivativesType1::dTdQ(data1D(indexQ), data1D(indexV), data1D(indexMIX), species_set, ND, dT(indexT));
		break;

	case 2:
		indexTe = 1;
		DerivativesType2::dTdQ(data1D(indexQ), data1D(indexV), data1D(indexMIX), species_set, ND, dT(indexT), dT(indexTe));
		break;

	case 3:
		indexTv = 1;
		DerivativesType3::dTdQ(data1D(indexQ), data1D(indexV), data1D(indexMIX), species_set, ND, dT(indexT), dT(indexTv));
		break;

	case 4:
		indexTv = 1;
		indexTe = 2;
		DerivativesType4::dTdQ(data1D(indexQ), data1D(indexV), data1D(indexMIX), species_set, ND, dT(indexT), dT(indexTv), dT(indexTe));
		break;

	case 5:
		indexTr = 1;
		DerivativesType5::dTdQ(data1D(indexQ), data1D(indexV), data1D(indexMIX), species_set, ND, dT(indexT), dT(indexTr));
		break;

	case 6:
		indexTr = 1;
		indexTe = 2;
		DerivativesType6::dTdQ(data1D(indexQ), data1D(indexV), data1D(indexMIX), species_set, ND, dT(indexT), dT(indexTr), dT(indexTe));
		break;

	case 7:
		indexTr = 1;
		indexTv = 2;
		DerivativesType7::dTdQ(data1D(indexQ), data1D(indexV), data1D(indexMIX), species_set, ND, dT(indexT), dT(indexTr), dT(indexTv));
		break;

	case 8:
		indexTr = 1;
		indexTv = 2;
		indexTe = 3;
		DerivativesType8::dTdQ(data1D(indexQ), data1D(indexV), data1D(indexMIX), species_set, ND, dT(indexT), dT(indexTr), dT(indexTv), dT(indexTe));
		break;
	}
}


void Derivatives::dpdQ(Data::DataStorageVector<Data::DataStorage>& data1D, Data::DataStorageVector<Data::DataStorage>& dT, CHEM::SpeciesSet& species_set, int ND,
				unsigned int type, unsigned int indexQ, unsigned int indexV, unsigned int indexW, unsigned int indexMIX,
				Data::DataStorage& dp)
{
	int indexT  = 0;
	int indexTe = 0;

	switch (type)
	{
	case 1:
		DerivativesType1::dpdQ(data1D(indexV), data1D(indexMIX), dT(indexT), species_set, ND, dp);
		break;

	case 2:
		indexTe = 1;
		DerivativesType2::dpdQ(data1D(indexV), data1D(indexMIX), dT(indexT), dT(indexTe), species_set, ND, dp);
		break;

	case 3:
		DerivativesType3::dpdQ(data1D(indexV), data1D(indexMIX), dT(indexT), species_set, ND, dp);
		break;

	case 4:
		indexTe = 2;
		DerivativesType4::dpdQ(data1D(indexV), data1D(indexMIX), dT(indexT), dT(indexTe), species_set, ND, dp);
		break;

	case 5:
		DerivativesType5::dpdQ(data1D(indexV), data1D(indexMIX), dT(indexT), species_set, ND, dp);
		break;

	case 6:
		indexTe = 2;
		DerivativesType6::dpdQ(data1D(indexV), data1D(indexMIX), dT(indexT), dT(indexTe), species_set, ND, dp);
		break;

	case 7:
		DerivativesType7::dpdQ(data1D(indexV), data1D(indexMIX), dT(indexT), species_set, ND, dp);
		break;

	case 8:
		indexTe = 3;
		DerivativesType8::dpdQ(data1D(indexV), data1D(indexMIX), dT(indexT), dT(indexTe), species_set, ND, dp);
		break;
	}
}


double Derivatives::a2(Data::DataStorageVector<Data::DataStorage>& data1D, Data::DataStorage& dp, CHEM::SpeciesSet& species_set, int ND,
				unsigned int type, unsigned int indexQ, unsigned int indexV, unsigned int indexW, unsigned int indexMIX)
{
	double o_a2;

	switch (type)
	{
	case 1:
		o_a2 = DerivativesType1::a2(data1D(indexQ), data1D(indexW), data1D(indexMIX), dp, species_set, ND);
		break;

	case 2:
		o_a2 = DerivativesType2::a2(data1D(indexQ), data1D(indexW), data1D(indexMIX), dp, species_set, ND);
		break;

	case 3:
		o_a2 = DerivativesType3::a2(data1D(indexQ), data1D(indexW), data1D(indexMIX), dp, species_set, ND);
		break;

	case 4:
		o_a2 = DerivativesType4::a2(data1D(indexQ), data1D(indexW), data1D(indexMIX), dp, species_set, ND);
		break;

	case 5:
		o_a2 = DerivativesType5::a2(data1D(indexQ), data1D(indexW), data1D(indexMIX), dp, species_set, ND);
		break;

	case 6:
		o_a2 = DerivativesType6::a2(data1D(indexQ), data1D(indexW), data1D(indexMIX), dp, species_set, ND);
		break;

	case 7:
		o_a2 = DerivativesType7::a2(data1D(indexQ), data1D(indexW), data1D(indexMIX), dp, species_set, ND);
		break;

	case 8:
		o_a2 = DerivativesType8::a2(data1D(indexQ), data1D(indexW), data1D(indexMIX), dp, species_set, ND);
		break;
	}

	return(o_a2);
}


}
}
