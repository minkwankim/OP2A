/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 22, 2015
 *      			Author: Minkwan Kim
 *
 * Variables.hpp
 * 			-  
 *  
 */
#ifndef VARIABLES_HPP_
#define VARIABLES_HPP_


#include "CFD/include/VariableChange.hpp"
#include "CFD/include/VariableConstants.hpp"


namespace OP2A{
namespace CFD{

class CFD_API CellData1D: public Common::NonInstantiable<CellData1D>
{
public:
	static void CompleteDataWIG(Data::DataStorageVector<Data::DataStorage>& cell1D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT,
								int indexQ, int indexV, int indexW, int indexMIX, int indexXs, int indexYs,
								int indexR, int indexdQ, int indexQnew, int indexS, int indexDivV, int indexdpdQ, int indexDs,int indexMus,
								int typeCase, bool is_initialize);

	static void CompleteDataWIG(Data::DataStorageVector<Data::DataStorage>& cell1D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, bool axis, bool viscous, int implicit, int typeCase, bool is_initialize);

	static void CompleteDataWIGCase1(Data::DataStorageVector<Data::DataStorage>& cell1D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, int indexQ, int indexV, int indexW, int indexMIX, int indexXs, int indexYs);
	static void CompleteDataWIGCase2(Data::DataStorageVector<Data::DataStorage>& cell1D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, int indexQ, int indexV, int indexW, int indexMIX, int indexXs, int indexYs);
	static void CompleteDataWIGCase3(Data::DataStorageVector<Data::DataStorage>& cell1D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, int indexQ, int indexV, int indexW, int indexMIX, int indexXs, int indexYs);
	static void CompleteDataWIGCase4(Data::DataStorageVector<Data::DataStorage>& cell1D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, int indexQ, int indexV, int indexW, int indexMIX, int indexXs, int indexYs);

	static void InitializeOtherData(Data::DataStorageVector<Data::DataStorage>& cell1D, unsigned int CFD_NT, int indexR, int indexdQ, int indexQnew, int indexS, int indexDivV, int indexdpdQ, int indexDs,int indexMus);
};


class CFD_API CellData2D: public Common::NonInstantiable<CellData2D>
{
public:
	static void InitializeData(Data::DataStorageVector<Data::DataStorage2D>& cell2D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, bool viscous, int implicit);
	static void InitializeData(Data::DataStorageVector<Data::DataStorage2D>& cell2D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, int indexdTdQ, int indexdSdQ, int indexKappas);

};


class CFD_API FaceData1D: public Common::NonInstantiable<FaceData1D>
{
public:
	static void InitializeData(Data::DataStorageVector<Data::DataStorage>& face1D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, bool viscous);

};

class CFD_API FaceData2D: public Common::NonInstantiable<FaceData2D>
{
public:
	static void InitializeData(Data::DataStorageVector<Data::DataStorage2D>& face2D, CHEM::SpeciesSet& species_set, int ND, unsigned int variabletype, unsigned int CFD_NT, bool viscous, int implicit);

};


}
}
#endif /* VARIABLES_HPP_ */
