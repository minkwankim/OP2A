/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 19, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChangeType2.cpp
 * 			-  
 *  
 */



#include <omp.h>

#include "CFD/include/VariableChangeTypes.hpp"


namespace OP2A{
namespace CFD{



void VariableChangeType2::V_to_Q(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_V(s);


	int indexT	= species_set.NS + ND;
	int indexTe	= indexT + 1;



	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_h0	+= data_V(s)*species_set.species[s].h0;
		}
	}


	// 2. E_k (Kinetic energy)
	double E_k = 0.0;
#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		E_k	+= rho_mix*data_V(k)*data_V(k);
	}
	E_k *= 0.5;


	// 3. E_int (translation/rotation energy) [Except electrons]
	double E_int = 0.0;
#pragma omp parallel for reduction(+:E_int)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_int	+= data_V(s) * (species_set.species[s].Cv_tr*data_V(indexT));
		}
	}

	// 4. E_e (Electron energy)
	double E_e = 0.0;
#pragma omp parallel for reduction(+:E_int)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			E_e	+= data_V(s) * (species_set.species[s].Cv_tra*data_V(indexTe) + species_set.species[s].h0);
		}
	}


	data_Q(indexT)	= E_int + E_k + E_h0 + E_e;
	data_Q(indexTe)	= E_e;
}



void VariableChangeType2::Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);


	int indexT	= species_set.NS + ND;
	int indexTe	= indexT + 1;


	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_h0	+= data_Q(s)*species_set.species[s].h0;
		}
	}


	// 2. E_k (Kinetic energy)
	double E_k = 0.0;
#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		E_k	+= data_Q(k)*data_Q(k) / rho_mix;
	}
	E_k *= 0.5;


	// 3. rhoCvmix
	double rhoCvmix = 0.0;
#pragma omp parallel for reduction(+:rhoCvmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCvmix	+= data_Q(s) * species_set.species[s].Cv_tr;
		}
	}

	// 4. rhoCve
	double rhoCve = 0.0;
#pragma omp parallel for reduction(+:rhoCve)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoCve	+= data_Q(s) * species_set.species[s].Cv_tra;
		}
	}

	double E_h0e = 0.0;
#pragma omp parallel for reduction(+:E_h0e)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			E_h0e	+= data_Q(s) * species_set.species[s].h0;
		}
	}




	double T, Te;
	T	= (data_Q(indexT) - data_Q(indexTe) - E_h0 - E_k) / rhoCvmix;
	Te	= (data_Q(indexTe) - E_h0e) / rhoCve;

	Common::ErrorCheckNonNegative<double>(T,  "Temperature(Type2)");
	Common::ErrorCheckNonNegative<double>(Te, "Electron Temperature(Type2)");

	data_V(indexT)	= T;
	data_V(indexTe)	= Te;
}




/*
 * V <=> W
 */
void VariableChangeType2::V_to_W(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTe	= indexT + 1;


	double rhoRmix	= 0.0;
#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoRmix	+= data_V(s) * species_set.species[s].R;
		}
	}


	double rhoeRe	= 0.0;
#pragma omp parallel for reduction(+:rhoeRe)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoeRe	+= data_V(s) * species_set.species[s].R;
		}
	}


	double pe	= rhoeRe * data_V(indexTe);
	double p	= rhoRmix * data_V(indexT) + pe;

	data_W(indexT)	= p;
	data_W(indexTe)	= pe;
}


void VariableChangeType2::W_to_V(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTe	= indexT + 1;


	double rhoRmix	= 0.0;
#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoRmix	+= data_W(s) * species_set.species[s].R;
		}
	}


	double rhoeRe	= 0.0;
#pragma omp parallel for reduction(+:rhoeRe)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoeRe	+= data_W(s) * species_set.species[s].R;
		}
	}


	double T	= (data_W(indexT) - data_W(indexTe)) / rhoRmix;
	double Te	= data_W(indexTe) / rhoeRe;

	Common::ErrorCheckNonNegative<double>(T,  "Temperature(Type2: W=>V)");
	data_V(indexT)	= T;
	data_V(indexTe)	= Te;
}


/*
 * Q <=> W
 */
void VariableChangeType2::Q_to_W(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTe	= indexT + 1;

	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);


	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_h0	+= data_Q(s)*species_set.species[s].h0;
		}
	}

	// 2. E_k (Kinetic energy)
	double E_k = 0.0;
#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		E_k	+= data_Q(k)*data_Q(k) / rho_mix;
	}
	E_k *= 0.5;

	// 3. rhoCvmix
	double rhoCvmix = 0.0;
#pragma omp parallel for reduction(+:rhoCvmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCvmix	+= data_Q(s) * species_set.species[s].Cv_tr;
		}
	}

	// 4. rhoCve
	double rhoCve = 0.0;
#pragma omp parallel for reduction(+:rhoCve)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoCve	+= data_Q(s) * species_set.species[s].Cv_tra;
		}
	}

	double E_h0e = 0.0;
#pragma omp parallel for reduction(+:E_h0e)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			E_h0e	+= data_Q(s) * species_set.species[s].h0;
		}
	}



	double T, Te;
	T	= (data_Q(indexT) - data_Q(indexTe) - E_h0 - E_k) / rhoCvmix;
	Te	= (data_Q(indexTe) - E_h0e) / rhoCve;

	Common::ErrorCheckNonNegative<double>(T,  "Temperature(Type2: Q =>W)");
	Common::ErrorCheckNonNegative<double>(Te, "Electron Temperature(Type2): Q => W");



	double rhoRmix	= 0.0;
#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoRmix	+= data_Q(s) * species_set.species[s].R;
		}
	}


	double rhoeRe	= 0.0;
#pragma omp parallel for reduction(+:rhoeRe)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoeRe	+= data_Q(s) * species_set.species[s].R;
		}
	}


	double pe	= rhoeRe * Te;
	double p	= rhoRmix * T + pe;

	data_W(indexT)	= p;
	data_W(indexTe)	= pe;
}

void VariableChangeType2::W_to_Q(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTe	= indexT + 1;


	double rhoRmix	= 0.0;
#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoRmix	+= data_W(s) * species_set.species[s].R;
		}
	}

	double rhoeRe	= 0.0;
#pragma omp parallel for reduction(+:rhoeRe)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoeRe	+= data_W(s) * species_set.species[s].R;
		}
	}


	double T	= (data_W(indexT) - data_W(indexTe)) / rhoRmix;
	double Te	= data_W(indexTe) / rhoeRe;

	Common::ErrorCheckNonNegative<double>(T,  "Temperature(Type2: W=>V)");




	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_W(s);



	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_h0	+= data_W(s)*species_set.species[s].h0;
		}
	}


	// 2. E_k (Kinetic energy)
	double E_k = 0.0;
#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		E_k	+= rho_mix*data_W(k)*data_W(k);
	}
	E_k *= 0.5;


	// 3. E_int (translation/rotation energy) [Except electrons]
	double E_int = 0.0;
#pragma omp parallel for reduction(+:E_int)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_int	+= data_W(s) * (species_set.species[s].Cv_tr*T);
		}
	}

	// 4. E_e (Electron energy)
	double E_e = 0.0;
#pragma omp parallel for reduction(+:E_int)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			E_e	+= data_W(s) * (species_set.species[s].Cv_tra*Te + species_set.species[s].h0);
		}
	}


	data_Q(indexT)	= E_int + E_h0 + E_k + E_e;
	data_Q(indexTe)	= E_e;
}





}
}
