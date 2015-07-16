/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 15, 2015
 *      			Author: Minkwan Kim
 *
 * Source_NONEQ.cpp
 * 			-  
 *  
 */

#include <omp.h>
#include <time.h>

#include "CFD/include/Source_NONEQ.hpp"
#include "Math/include/MathMisc.hpp"

namespace OP2A{
namespace CFD{

void Source_NONEQ::S_NONEQ(Data::DataStorageVector<Data::DataStorage>& data1D, CHEM::SpeciesSet& species_set, int ND,
		unsigned int type, unsigned int indexQ, unsigned int indexV, unsigned int indexW, unsigned int indexS)
{
	int indexTr;
	int indexTv;
	int indexTe;
	int indexT	= species_set.NS + ND;

	switch (type)
	{
	case 2:
		indexTe = indexT + 1;
		data1D(indexS)(indexTe)	= S_he(data1D, species_set, ND, indexTe, indexQ, indexV, indexW);
		break;
	}
}




}
}
