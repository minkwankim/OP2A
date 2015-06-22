/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 16, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChangeMixture.hpp
 * 			-  
 *  
 */


#ifndef VARIABLECHANGEMIXTURE_HPP_
#define VARIABLECHANGEMIXTURE_HPP_

#include <limits>
#include "CFD/include/CFDCommon.hpp"


namespace OP2A{
namespace CFD{


class CFD_API VariableChangeMixture : public Common::NonInstantiable<VariableChangeMixture>
{
public:
	static void Xs(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Xs, unsigned int CFD_NT);
	static void Ys(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Ys, unsigned int CFD_NT);

	static double rho_mix(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, unsigned int CFD_NT);
	static double R_mix(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, unsigned int CFD_NT);
	static double M_mix(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, unsigned int CFD_NT);
	static double Cv_tra_mix(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, unsigned int CFD_NT);
	static double Cv_rot_mix(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, unsigned int CFD_NT);
	static double Cv_tr_mix(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, unsigned int CFD_NT);

	static void MIX(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Ys,  Data::DataStorage& data_Xs, unsigned int type, unsigned int CFD_NT, Data::DataStorage& data_MIX);
};


}
}



#endif /* VARIABLECHANGEMIXTURE_HPP_ */
