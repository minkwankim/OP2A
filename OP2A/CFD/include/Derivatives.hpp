/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 25, 2015
 *      			Author: Minkwan Kim
 *
 * dTdQ.hpp
 * 			-  
 *  
 */
#ifndef DTDQ_HPP_
#define DTDQ_HPP_


#include "CFD/include/CFDCommon.hpp"
#include "DATA/include/DataStorage2D.hpp"
#include "DATA/include/DataStorageVector.hpp"



namespace OP2A{
namespace CFD{

class CFD_API DerivativesType1 : public Common::NonInstantiable<DerivativesType1>
{
public:
	static void dTdQ(Data::DataStorage& data_Q, Data::DataStorage& data_V, 		Data::DataStorage& data_MIX,	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT);
	static void dpdQ(Data::DataStorage& data_V, Data::DataStorage& data_MIX,	Data::DataStorage& dT, 			CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dp);
	static double a2(Data::DataStorage& data_Q, Data::DataStorage& data_W,  	Data::DataStorage& data_MIX, 	Data::DataStorage& dp, CHEM::SpeciesSet& species_set, int ND);

	static void d2TdQ2(Data::DataStorage& data_Q, Data::DataStorage& data_V,  	Data::DataStorage& data_MIX,	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT, Data::DataStorage2D& dT2);
};

class CFD_API DerivativesType2 : public Common::NonInstantiable<DerivativesType2>
{
public:
	static void dTdQ(Data::DataStorage& data_Q, Data::DataStorage& data_V,		Data::DataStorage& data_MIX,	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT,  Data::DataStorage& dTe);
	static void dpdQ(Data::DataStorage& data_V, Data::DataStorage& data_MIX,	Data::DataStorage& dT, 			Data::DataStorage& dTe, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dp);
	static double a2(Data::DataStorage& data_Q, Data::DataStorage& data_W, 		Data::DataStorage& data_MIX,	Data::DataStorage& dp, CHEM::SpeciesSet& species_set, int ND);

	static void d2TdQ2(Data::DataStorage& data_Q, Data::DataStorage& data_V,  	Data::DataStorage& data_MIX,	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT, Data::DataStorage& dTe, Data::DataStorage2D& d2T, Data::DataStorage2D& d2Te);
};


class CFD_API DerivativesType3 : public Common::NonInstantiable<DerivativesType3>
{
public:
	static void dTdQ(Data::DataStorage& data_Q, Data::DataStorage& data_V, 		Data::DataStorage& data_MIX, 	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT,  Data::DataStorage& dTv);
	static void dpdQ(Data::DataStorage& data_V, Data::DataStorage& data_MIX,	Data::DataStorage& dT, 			CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dp);
	static double a2(Data::DataStorage& data_Q, Data::DataStorage& data_W, 		Data::DataStorage& data_MIX,	Data::DataStorage& dp, CHEM::SpeciesSet& species_set, int ND);

	static void d2TdQ2(Data::DataStorage& data_Q, Data::DataStorage& data_V,  	Data::DataStorage& data_MIX,	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT, Data::DataStorage2D& dT2);
};

class CFD_API DerivativesType4 : public Common::NonInstantiable<DerivativesType4>
{
public:
	static void dTdQ(Data::DataStorage& data_Q, Data::DataStorage& data_V, 		Data::DataStorage& data_MIX,	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT, Data::DataStorage& dTv, Data::DataStorage& dTe);
	static void dpdQ(Data::DataStorage& data_V, Data::DataStorage& data_MIX,	Data::DataStorage& dT, 			Data::DataStorage& dTe, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dp);
	static double a2(Data::DataStorage& data_Q, Data::DataStorage& data_W, 		Data::DataStorage& data_MIX,	Data::DataStorage& dp, CHEM::SpeciesSet& species_set, int ND);

	static void d2TdQ2(Data::DataStorage& data_Q, Data::DataStorage& data_V,  	Data::DataStorage& data_MIX,	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT, Data::DataStorage& dTe, Data::DataStorage2D& d2T, Data::DataStorage2D& d2Te);
};

class CFD_API DerivativesType5 : public Common::NonInstantiable<DerivativesType5>
{
public:
	static void dTdQ(Data::DataStorage& data_Q, Data::DataStorage& data_V, 		Data::DataStorage& data_MIX, 	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT,  Data::DataStorage& dTr);
	static void dpdQ(Data::DataStorage& data_V, Data::DataStorage& data_MIX, 	Data::DataStorage& dT, 			CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dp);
	static double a2(Data::DataStorage& data_Q, Data::DataStorage& data_W, 		Data::DataStorage& data_MIX, 	Data::DataStorage& dp, CHEM::SpeciesSet& species_set, int ND);

	static void d2TdQ2(Data::DataStorage& data_Q, Data::DataStorage& data_V,  	Data::DataStorage& data_MIX,	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT, Data::DataStorage2D& dT2);
};

class CFD_API DerivativesType6 : public Common::NonInstantiable<DerivativesType6>
{
public:
	static void dTdQ(Data::DataStorage& data_Q, Data::DataStorage& data_V, 		Data::DataStorage& data_MIX,	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT,  Data::DataStorage& dTr,   Data::DataStorage& dTe);
	static void dpdQ(Data::DataStorage& data_V, Data::DataStorage& data_MIX,	Data::DataStorage& dT,			Data::DataStorage& dTe, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dp);
	static double a2(Data::DataStorage& data_Q, Data::DataStorage& data_W, 		Data::DataStorage& data_MIX,	Data::DataStorage& dp, CHEM::SpeciesSet& species_set, int ND);

	static void d2TdQ2(Data::DataStorage& data_Q, Data::DataStorage& data_V,  	Data::DataStorage& data_MIX,	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT, Data::DataStorage& dTe, Data::DataStorage2D& d2T, Data::DataStorage2D& d2Te);
};

class CFD_API DerivativesType7 : public Common::NonInstantiable<DerivativesType7>
{
public:
	static void dTdQ(Data::DataStorage& data_Q, Data::DataStorage& data_V, 		Data::DataStorage& data_MIX, 	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT,  Data::DataStorage& dTr,   Data::DataStorage& dTv);
	static void dpdQ(Data::DataStorage& data_V, Data::DataStorage& data_MIX,	Data::DataStorage& dT, 			CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dp);
	static double a2(Data::DataStorage& data_Q, Data::DataStorage& data_W, 		Data::DataStorage& data_MIX,	Data::DataStorage& dp, CHEM::SpeciesSet& species_set, int ND);

	static void d2TdQ2(Data::DataStorage& data_Q, Data::DataStorage& data_V,  	Data::DataStorage& data_MIX,	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT, Data::DataStorage2D& dT2);
};

class CFD_API DerivativesType8 : public Common::NonInstantiable<DerivativesType8>
{
public:
	static void dTdQ(Data::DataStorage& data_Q, Data::DataStorage& data_V, 		Data::DataStorage& data_MIX,	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT,  Data::DataStorage& dTr,   Data::DataStorage& dTv, Data::DataStorage& dTe);
	static void dpdQ(Data::DataStorage& data_V, Data::DataStorage& data_MIX,	Data::DataStorage& dT, 			Data::DataStorage& dTe, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dp);
	static double a2(Data::DataStorage& data_Q, Data::DataStorage& data_W, 		Data::DataStorage& data_MIX,	Data::DataStorage& dp, CHEM::SpeciesSet& species_set, int ND);

	static void d2TdQ2(Data::DataStorage& data_Q, Data::DataStorage& data_V,  	Data::DataStorage& data_MIX,	CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& dT, Data::DataStorage& dTe, Data::DataStorage2D& d2T, Data::DataStorage2D& d2Te);
};



class CFD_API Derivatives : public Common::NonInstantiable<Derivatives>
{
public:
	static void dTdQ(Data::DataStorageVector<Data::DataStorage>& data1D, CHEM::SpeciesSet& species_set, int ND,
					unsigned int type, unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX,
					Data::DataStorageVector<Data::DataStorage>& dT);

	static void dpdQ(Data::DataStorageVector<Data::DataStorage>& data1D, Data::DataStorageVector<Data::DataStorage>& dT, CHEM::SpeciesSet& species_set, int ND,
					unsigned int type, unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX,
					Data::DataStorage& dp);

	static double a2(Data::DataStorageVector<Data::DataStorage>& data1D, Data::DataStorage& dp, CHEM::SpeciesSet& species_set, int ND,
					unsigned int type, unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX);

	static void d2pdQ(Data::DataStorageVector<Data::DataStorage>& data1D,
					  Data::DataStorageVector<Data::DataStorage>& dT,
					  CHEM::SpeciesSet& species_set, int ND,
					  unsigned int type, unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX,
					  Data::DataStorage2D& d2p);

	static void d2pdQ_CaseA1(Data::DataStorageVector<Data::DataStorage>& data1D,
							Data::DataStorage& dT,
							CHEM::SpeciesSet& species_set, int ND,
							unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX,
							Data::DataStorage2D& d2p);

	static void d2pdQ_CaseA2(Data::DataStorageVector<Data::DataStorage>& data1D,
								Data::DataStorage& dT,
								CHEM::SpeciesSet& species_set, int ND,
								unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX,
								Data::DataStorage2D& d2p);

	static void d2pdQ_CaseA(Data::DataStorageVector<Data::DataStorage>& data1D,
							Data::DataStorage& dT,
							Data::DataStorage2D& d2T,
							CHEM::SpeciesSet& species_set, int ND,
							unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX,
							Data::DataStorage2D& d2p);

	static void d2pdQ_CaseB(Data::DataStorageVector<Data::DataStorage>& data1D,
							Data::DataStorage& dT,
							Data::DataStorage& dTe,
							Data::DataStorage2D& d2T,
							Data::DataStorage2D& d2Te,
							CHEM::SpeciesSet& species_set, int ND,
							unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX,
							Data::DataStorage2D& d2p);

	static void da2dQ(Data::DataStorageVector<Data::DataStorage>& data1D,
						Data::DataStorage& dp,
						Data::DataStorage2D& d2p,
						CHEM::SpeciesSet& species_set, int ND,
						unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX,
						Data::DataStorage& da2);

	static void dadQ(Data::DataStorage& da2, double a, Data::DataStorage& da);

};



}
}


#endif /* DTDQ_HPP_ */
