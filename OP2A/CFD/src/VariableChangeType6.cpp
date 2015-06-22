/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 19, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChangeType6.cpp
 * 			-  
 *  
 */

#include <omp.h>

#include "CFD/include/VariableChangeTypes.hpp"


namespace OP2A{
namespace CFD{



void VariableChangeType6::V_to_Q(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;
	int indexTe	= species_set.NS + ND + 2;

	double T	= data_V(indexT);
	double Tr	= data_V(indexTr);
	double Te	= data_V(indexTe);

	double E, Er, Ee;



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


	// 3. E_tra (translation energy) [Except electrons]
	double E_tra = 0.0;
#pragma omp parallel for reduction(+:E_tra)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_tra	+= data_V(s) * (species_set.species[s].Cv_tra*T);
		}
	}


	// 4. Erot (translation energy) [Except electrons]
	Er = 0.0;
#pragma omp parallel for reduction(+:Er)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			Er	+= data_V(s) * (species_set.species[s].Cv_rot*Tr);
		}
	}


	// 5. E_e (Electron energy)
	Ee = 0.0;
#pragma omp parallel for reduction(+:Ee)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			Ee	+= data_V(s) * (species_set.species[s].Cv_tra*Te + species_set.species[s].h0);
		}
	}


	E	= E_tra + E_h0 + E_k + Er + Ee;

	data_Q(indexT)	= E;
	data_Q(indexTr)	= Er;
	data_Q(indexTe)	= Ee;
}

void VariableChangeType6::Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;
	int indexTe	= species_set.NS + ND + 2;

	double E	= data_Q(indexT);
	double Er	= data_Q(indexTr);
	double Ee	= data_Q(indexTe);

	double T, Tr, Te;



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


	// 3. rhoCv_tra_mix
	double rhoCv_tra_mix = 0.0;
#pragma omp parallel for reduction(+:rhoCv_tra_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCv_tra_mix	+= data_Q(s) * species_set.species[s].Cv_tra;
		}
	}


	// 4. rhoCv_rot_mix
	double rhoCv_rot_mix = 0.0;
#pragma omp parallel for reduction(+:rhoCv_rot_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCv_rot_mix	+= data_Q(s) * species_set.species[s].Cv_rot;
		}
	}


	// 5. rhoCve
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


	T	= (E - Er - Ee - E_h0 - E_k) / rhoCv_tra_mix;
	Tr	= Er / rhoCv_rot_mix;
	Te	= (Ee - E_h0e) / rhoCve;

	Common::ErrorCheckNonNegative<double>(T,  "Temperature(Type6): Q => V");
	Common::ErrorCheckNonNegative<double>(Tr,  "Rotational Temperature(Type6): Q => V");
	Common::ErrorCheckNonNegative<double>(Te, "Electron Temperature(Type6): Q => V");

	data_V(indexT)	= T;
	data_V(indexTr)	= Tr;
	data_V(indexTe)	= Te;
}



/*
 * V <=> W
 */
void VariableChangeType6::V_to_W(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;
	int indexTe	= species_set.NS + ND + 2;

	double T	= data_V(indexT);
	double Tr	= data_V(indexTr);
	double Te	= data_V(indexTe);

	double p, pe;



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


	pe	= rhoeRe * Te;
	p	= rhoRmix * T + pe;

	data_W(indexT)	= p;
	data_W(indexTr)	= Tr;
	data_W(indexTe)	= pe;
}


void VariableChangeType6::W_to_V(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;
	int indexTe	= species_set.NS + ND + 2;

	double p	= data_W(indexT);
	double Tr	= data_W(indexTr);
	double pe	= data_W(indexTe);

	double T, Te;


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


	T	= (p - pe) / rhoRmix;
	Te	= pe / rhoeRe;

	Common::ErrorCheckNonNegative<double>(T,  "Temperature(Type6: W => V)");

	data_V(indexT)	= T;
	data_V(indexTr)	= Tr;
	data_V(indexTe)	= Te;
}




/*
 * Q <=> W
 */
void VariableChangeType6::Q_to_W(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;
	int indexTe	= species_set.NS + ND + 2;

	double E	= data_Q(indexT);
	double Er	= data_Q(indexTr);
	double Ee	= data_Q(indexTe);

	double T, Tr, Te;



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


	// 3. rhoCv_tra_mix
	double rhoCv_tra_mix = 0.0;
#pragma omp parallel for reduction(+:rhoCv_tra_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCv_tra_mix	+= data_Q(s) * species_set.species[s].Cv_tra;
		}
	}


	// 4. rhoCv_rot_mix
	double rhoCv_rot_mix = 0.0;
#pragma omp parallel for reduction(+:rhoCv_rot_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCv_rot_mix	+= data_Q(s) * species_set.species[s].Cv_rot;
		}
	}


	// 5. rhoCve
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


	T	= (E - Er - Ee - E_h0 - E_k) / rhoCv_tra_mix;
	Tr	= Er / rhoCv_rot_mix;
	Te	= (Ee - E_h0e) / rhoCve;

	Common::ErrorCheckNonNegative<double>(T,  "Temperature(Type6): Q => V");
	Common::ErrorCheckNonNegative<double>(Tr,  "Rotational Temperature(Type6): Q => V");
	Common::ErrorCheckNonNegative<double>(Te, "Electron Temperature(Type6): Q => V");



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
	data_W(indexTr)	= Tr;
	data_W(indexTe)	= pe;
}






void VariableChangeType6::W_to_Q(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;
	int indexTe	= species_set.NS + ND + 2;

	double p	= data_W(indexT);
	double Tr	= data_W(indexTr);
	double pe	= data_W(indexTe);

	double T, Te;


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


	T	= (p - pe) / rhoRmix;
	Te	= pe / rhoeRe;

	Common::ErrorCheckNonNegative<double>(T,  "Temperature(Type6: W => V)");


	double E, Er, Ee;

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


	// 3. E_tra (translation energy) [Except electrons]
	double E_tra = 0.0;
#pragma omp parallel for reduction(+:E_tra)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_tra	+= data_W(s) * (species_set.species[s].Cv_tra*T);
		}
	}


	// 4. Erot (translation energy) [Except electrons]
	Er = 0.0;
#pragma omp parallel for reduction(+:Er)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			Er	+= data_W(s) * (species_set.species[s].Cv_rot*Tr);
		}
	}


	// 5. E_e (Electron energy)
	Ee = 0.0;
#pragma omp parallel for reduction(+:Ee)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			Ee	+= data_W(s) * (species_set.species[s].Cv_tra*Te + species_set.species[s].h0);
		}
	}


	E	= E_tra + E_h0 + E_k + Er + Ee;


	data_Q(indexT)	= E;
	data_Q(indexTr)	= Er;
	data_Q(indexTe)	= Ee;
}












}
}
