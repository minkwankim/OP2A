/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 16, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChangeMixture.cpp
 * 			-  
 *  
 */



#include <omp.h>

#include "CFD/include/VariableChangeMixture.hpp"
#include "CFD/include/VariableSet.hpp"



namespace OP2A{
namespace CFD{

// Mole fraction
void VariableChangeMixture::Xs(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Xs, unsigned int CFD_NT)
{
#pragma omp parallel for num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	data_Xs(s)	= data_Q(s) / species_set.species[s].m;


	// Total number of Moles
	double n_mix	= 0.0;
#pragma omp parallel for reduction(+:n_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	n_mix	+= data_Xs(s);

#pragma omp parallel for num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	data_Xs(s)	/= n_mix;
}


// Mass fraction
void VariableChangeMixture::Ys(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Ys, unsigned int CFD_NT)
{
	// Total number of Moles
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);

#pragma omp parallel for num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	data_Ys(s) = data_Q(s)/rho_mix;
}

/*
 * Mixture properties
 * @author Minkwan Kim
 * @version 22/6/2015 1.0
 */
double VariableChangeMixture::rho_mix(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, unsigned int CFD_NT)
{
	double o_rhomix = 0.0;
#pragma omp parallel for reduction(+:o_rhomix)
	for (int s = 0; s <= species_set.NS-1; s++)	o_rhomix	+= data_Q(s);

	return o_rhomix;
}


double VariableChangeMixture::R_mix(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, unsigned int CFD_NT)
{
	double o_Rmix = 0.0;
#pragma omp parallel for reduction(+:o_Rmix)
	for (int s = 0; s <= species_set.NS-1; s++)	o_Rmix	+= data_Ys(s) * species_set.species[s].R;

	return o_Rmix;
}

double VariableChangeMixture::M_mix(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, unsigned int CFD_NT)
{
	double o_Mmix = 0.0;
#pragma omp parallel for reduction(+:o_Mmix)
	for (int s = 0; s <= species_set.NS-1; s++)	o_Mmix	+= data_Xs(s) * species_set.species[s].M;

	return o_Mmix;
}



double VariableChangeMixture::Cv_tra_mix(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, unsigned int CFD_NT)
{
	double o_Cv_tra_mix = 0.0;
#pragma omp parallel for reduction(+:o_Cv_tra_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	o_Cv_tra_mix	+= data_Ys(s) * species_set.species[s].Cv_tra;

	return o_Cv_tra_mix;
}

double VariableChangeMixture::Cv_rot_mix(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, unsigned int CFD_NT)
{
	double o_Cv_tra_mix = 0.0;
#pragma omp parallel for reduction(+:o_Cv_tra_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	o_Cv_tra_mix	+= data_Ys(s) * species_set.species[s].Cv_rot;

	return o_Cv_tra_mix;
}

double VariableChangeMixture::Cv_tr_mix(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, unsigned int CFD_NT)
{
	double o_Cv_tra_mix = 0.0;
#pragma omp parallel for reduction(+:o_Cv_tra_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	o_Cv_tra_mix	+= data_Ys(s) * species_set.species[s].Cv_tr;

	return o_Cv_tra_mix;
}




void VariableChangeMixture::MIX(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, unsigned int type, unsigned int CFD_NT, Data::DataStorage& data_Ys,  Data::DataStorage& data_Xs, Data::DataStorage& data_MIX)
{
	Xs(data_Q, species_set, ND, data_Xs, CFD_NT);
	Ys(data_Q, species_set, ND, data_Ys, CFD_NT);


	int index;
	double temp_rhoMix;
	double temp_RMix;
	double temp_MMix;
	double temp_CvtraMix;
	double temp_CvrotMix;
	double temp_CvtrMix;


	temp_rhoMix	= rho_mix(data_Q, species_set, ND, data_Xs, data_Ys, CFD_NT);
	temp_RMix	= R_mix(data_Q, species_set, ND, data_Xs, data_Ys, CFD_NT);
	temp_MMix	= M_mix(data_Q, species_set, ND, data_Xs, data_Ys, CFD_NT);

	if (type <= 4)
	{
		temp_CvtrMix	= Cv_tr_mix(data_Q, species_set, ND, data_Xs, data_Ys, CFD_NT);

		index	= data_MIX.dataMap.find(CFD_VariableSet::rho_mix());
		data_MIX(index)	= temp_rhoMix;

		index	= data_MIX.dataMap.find(CFD_VariableSet::R_mix());
		data_MIX(index)	= temp_RMix;

		index	= data_MIX.dataMap.find(CFD_VariableSet::M_mix());
		data_MIX(index)	= temp_MMix;

		index	= data_MIX.dataMap.find(CFD_VariableSet::Cv_mix(CHEM::EnergyMode::E_TR));
		data_MIX(index)	= temp_CvtrMix;

		index	= data_MIX.dataMap.find(CFD_VariableSet::gamma_mix());
		data_MIX(index)	= (temp_CvtrMix + temp_RMix) / temp_CvtrMix;
	}
	else
	{
		temp_CvtraMix	= Cv_tra_mix(data_Q, species_set, ND, data_Xs, data_Ys, CFD_NT);
		temp_CvrotMix	= Cv_rot_mix(data_Q, species_set, ND, data_Xs, data_Ys, CFD_NT);

		index	= data_MIX.dataMap.find(CFD_VariableSet::rho_mix());
		data_MIX(index)	= temp_rhoMix;

		index	= data_MIX.dataMap.find(CFD_VariableSet::R_mix());
		data_MIX(index)	= temp_RMix;

		index	= data_MIX.dataMap.find(CFD_VariableSet::M_mix());
		data_MIX(index)	= temp_MMix;

		index	= data_MIX.dataMap.find(CFD_VariableSet::Cv_mix(CHEM::EnergyMode::E_TRA));
		data_MIX(index)	= temp_CvtraMix;

		index	= data_MIX.dataMap.find(CFD_VariableSet::Cv_mix(CHEM::EnergyMode::E_ROT));
		data_MIX(index)	= temp_CvrotMix;


		index	= data_MIX.dataMap.find(CFD_VariableSet::gamma_mix());
		data_MIX(index)	= (temp_CvtraMix + temp_CvrotMix + temp_RMix) / (temp_CvtraMix + temp_CvrotMix);
	}
}





}
}
