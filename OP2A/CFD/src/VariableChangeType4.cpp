/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 19, 2015
 *      			Author: Minkwan Kim
 *
 * VariavleChangeType4.cpp
 * 			-  
 *  
 */



#include <omp.h>

#include "CFD/include/VariableChangeTypes.hpp"


namespace OP2A{
namespace CFD{



void VariableChangeType4::V_to_Q(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTve	= indexT + 1;
	int indexTe		= indexT + 2;

	double T	= data_V(indexT);
	double Tve	= data_V(indexTve);
	double Te	= data_V(indexTe);

	double E;
	double Eve;
	double Ee;



	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_h0	+= data_V(s)*species_set.species[s].h0;
		}
	}


	// 2. E_k (Kinetic energy)
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_V(s);

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
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_int	+= data_V(s) * (species_set.species[s].Cv_tr*T);
		}
	}


	// 4. E_VE (VIBRATION-ELECTRONIC energy)
	Eve = 0.0;
#pragma omp parallel for reduction(+:Eve)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			Eve	+= data_V(s) * species_set.species[s].e_VE(Tve);
		}
	}


	// 5. Ee
	Ee = 0.0;
#pragma omp parallel for reduction(+:Ee)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			Ee	+= data_V(s) * (species_set.species[s].Cv_tra*Te + species_set.species[s].h0);
		}
	}

	E	= E_int + E_h0 + E_k + Eve + Ee;

	data_Q(indexT)		= E;
	data_Q(indexTve)	= Eve;
	data_Q(indexTe)		= Ee;
}

void VariableChangeType4::Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTve	= indexT + 1;
	int indexTe		= indexT + 2;

	double E	= data_Q(indexT);
	double Eve	= data_Q(indexTve);
	double Ee	= data_Q(indexTe);

	double T;
	double Tve;
	double Te;



	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_h0	+= data_Q(s)*species_set.species[s].h0;
		}
	}


	// 2. E_k (Kinetic energy)
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);

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
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCvmix	+= data_Q(s) * species_set.species[s].Cv_tr;
		}
	}


	// 4. rhoCve
	double rhoCve = 0.0;
#pragma omp parallel for reduction(+:rhoCve)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoCve	+= data_Q(s) * species_set.species[s].Cv_tra;
		}
	}

	double E_h0e = 0.0;
#pragma omp parallel for reduction(+:E_h0e)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			E_h0e	+= data_Q(s) * species_set.species[s].h0;
		}
	}



	T	= (E - E_h0 - E_k - Eve - Ee) / rhoCvmix;
	Tve	= CFD_calculate_Tve(data_Q, species_set, Eve, 10000, 0.0001, CFD_NT);
	Te	= (Ee - E_h0e) / rhoCve;

	Common::ErrorCheckNonNegative<double>(T,   "Temperature(Type4): Q => V");
	Common::ErrorCheckNonNegative<double>(Tve, "Vibrational Temperature(Type4): Q => V");
	Common::ErrorCheckNonNegative<double>(Te,  "Electron Temperature(Type4): Q => V");

	data_V(indexT)		= T;
	data_V(indexTve)	= Tve;
	data_V(indexTe)	= Te;
}




//////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * V ><=> W
 */
void VariableChangeType4::V_to_W(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTve	= species_set.NS + ND + 1;
	int indexTe		= species_set.NS + ND + 2;

	double T	= data_V(indexT);
	double Tve	= data_V(indexTve);
	double Te	= data_V(indexTe);



	double rhoRmix	= 0.0;
#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoRmix	+= data_V(s)*species_set.species[s].R;
		}
	}

	double rhoeRe	= 0.0;
#pragma omp parallel for reduction(+:rhoeRe)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoeRe	+= data_V(s)*species_set.species[s].R;
		}
	}


	double pe	= rhoeRe * Te;
	double p	= rhoRmix*T + pe;

	data_W(indexT)		= p;
	data_W(indexTve)	= Tve;
	data_W(indexTe)	= pe;
}

void VariableChangeType4::W_to_V(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTve	= species_set.NS + ND + 1;
	int indexTe		= species_set.NS + ND + 2;

	double p	= data_W(indexT);
	double Tve	= data_W(indexTve);
	double pe	= data_W(indexTe);



	double rhoRmix	= 0.0;
#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoRmix	+= data_W(s)*species_set.species[s].R;
		}
	}

	double rhoeRe	= 0.0;
#pragma omp parallel for reduction(+:rhoeRe)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoeRe	+= data_W(s)*species_set.species[s].R;
		}
	}

	double Te	= pe / rhoeRe;
	double T	= (p - pe) / rhoRmix;

	data_V(indexT)		= T;
	data_V(indexTve)	= Tve;
	data_V(indexTe)		= Te;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Q <=> W
 */
void VariableChangeType4::Q_to_W(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTve	= indexT + 1;
	int indexTe		= indexT + 2;

	double E	= data_Q(indexT);
	double Eve	= data_Q(indexTve);
	double Ee	= data_Q(indexTe);

	double T;
	double Tve;
	double Te;
	double p;
	double pe;




	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_h0	+= data_Q(s)*species_set.species[s].h0;
		}
	}


	// 2. E_k (Kinetic energy)
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);

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
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCvmix	+= data_Q(s) * species_set.species[s].Cv_tr;
		}
	}


	// 4. rhoCve
	double rhoCve = 0.0;
#pragma omp parallel for reduction(+:rhoCve)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoCve	+= data_Q(s) * species_set.species[s].Cv_tra;
		}
	}

	double E_h0e = 0.0;
#pragma omp parallel for reduction(+:E_h0e)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			E_h0e	+= data_Q(s) * species_set.species[s].h0;
		}
	}



	T	= (E - E_h0 - E_k - Eve - Ee) / rhoCvmix;
	Tve	= CFD_calculate_Tve(data_Q, species_set, Eve, 10000, 0.0001, CFD_NT);
	Te	= (Ee - E_h0e) / rhoCve;

	Common::ErrorCheckNonNegative<double>(T,   "Temperature(Type4): Q => V");
	Common::ErrorCheckNonNegative<double>(Tve, "Vibrational Temperature(Type4): Q => V");
	Common::ErrorCheckNonNegative<double>(Te,  "Electron Temperature(Type4): Q => V");


	double rhoRmix	= 0.0;
#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoRmix	+= data_Q(s)*species_set.species[s].R;
		}
	}

	double rhoeRe	= 0.0;
#pragma omp parallel for reduction(+:rhoeRe)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoeRe	+= data_Q(s)*species_set.species[s].R;
		}
	}


	pe	= rhoeRe * Te;
	p	= rhoRmix*T + pe;

	data_W(indexT)		= p;
	data_W(indexTve)	= Tve;
	data_W(indexTe)	= pe;
}

void VariableChangeType4::W_to_Q(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTve	= species_set.NS + ND + 1;
	int indexTe		= species_set.NS + ND + 2;

	double p	= data_W(indexT);
	double pe	= data_W(indexTe);


	double rhoRmix	= 0.0;
#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoRmix	+= data_W(s)*species_set.species[s].R;
		}
	}

	double rhoeRe	= 0.0;
#pragma omp parallel for reduction(+:rhoeRe)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoeRe	+= data_W(s)*species_set.species[s].R;
		}
	}

	double Te	= pe / rhoeRe;
	double T	= (p - pe) / rhoRmix;
	double Tve	= data_W(indexTve);


	double E;
	double Eve;
	double Ee;



	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_h0	+= data_W(s)*species_set.species[s].h0;
		}
	}


	// 2. E_k (Kinetic energy)
	double rho_mix	= 0.0;
#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_W(s);

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
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_int	+= data_W(s) * (species_set.species[s].Cv_tr*T);
		}
	}


	// 4. E_VE (VIBRATION-ELECTRONIC energy)
	Eve = 0.0;
#pragma omp parallel for reduction(+:Eve)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			Eve	+= data_W(s) * species_set.species[s].e_VE(Tve);
		}
	}


	// 5. Ee
	Ee = 0.0;
#pragma omp parallel for reduction(+:Ee)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			Ee	+= data_W(s) * (species_set.species[s].Cv_tra*Te + species_set.species[s].h0);
		}
	}

	E	= E_int + E_h0 + E_k + Eve + Ee;

	data_Q(indexT)		= E;
	data_Q(indexTve)	= Eve;
	data_Q(indexTe)		= Ee;
}

}
}
