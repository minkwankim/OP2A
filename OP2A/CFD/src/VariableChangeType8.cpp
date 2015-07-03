/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 19, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChangeType8.cpp
 * 			-  
 *  
 */



#include <omp.h>

#include "CFD/include/VariableChangeTypes.hpp"


namespace OP2A{
namespace CFD{


void VariableChangeType8::V_to_Q(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTr		= indexT + 1;
	int indexTve	= indexT + 2;
	int indexTe		= indexT + 3;

	double T	= data_V(indexT);
	double Tr	= data_V(indexTr);
	double Tve	= data_V(indexTve);
	double Te	= data_V(indexTe);

	double E;
	double Er;
	double Eve;
	double Ee;



	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_h0	+= data_V(s)*species_set.species[s].h0;
		}
	}


	// 2. E_k (Kinetic energy)
	double rho_mix	= 0.0;
//#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_V(s);

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		E_k	+= rho_mix*data_V(k)*data_V(k);
	}
	E_k *= 0.5;


	// 3. E_tra (translation/rotation energy) [Except electrons]
	double E_tra = 0.0;
//#pragma omp parallel for reduction(+:E_tra)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_tra	+= data_V(s) * (species_set.species[s].Cv_tra*T);
		}
	}


	// 4. E_rot (translation/rotation energy) [Except electrons]
	Er = 0.0;
//#pragma omp parallel for reduction(+:Er)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			Er	+= data_V(s) * (species_set.species[s].Cv_rot*Tr);
		}
	}


	// 5. E_VE (VIBRATION-ELECTRONIC energy)
	Eve = 0.0;
//#pragma omp parallel for reduction(+:Eve)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			Eve	+= data_V(s) * species_set.species[s].e_VE(Tve);
		}
	}


	// 5. Ee
	Ee = 0.0;
//#pragma omp parallel for reduction(+:Ee)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			Ee	= data_V(s) * (species_set.species[s].Cv_tra*Te + species_set.species[s].h0);
			break;
		}
	}

	E	= E_tra + E_h0 + E_k + Er + Eve + Ee;

	data_Q(indexT)		= E;
	data_Q(indexTr)		= Er;
	data_Q(indexTve)	= Eve;
	data_Q(indexTe)		= Ee;
}

void VariableChangeType8::Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTr	= indexT + 1;
	int indexTve	= indexT + 2;
	int indexTe		= indexT + 3;

	double E	= data_Q(indexT);
	double Er	= data_Q(indexTr);
	double Eve	= data_Q(indexTve);
	double Ee	= data_Q(indexTe);

	double T;
	double Tr;
	double Tve;
	double Te;



	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_h0	+= data_Q(s)*species_set.species[s].h0;
		}
	}


	// 2. E_k (Kinetic energy)
	double rho_mix	= 0.0;
//#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		E_k	+= data_Q(k)*data_Q(k) / rho_mix;
	}
	E_k *= 0.5;


	// 3. rhoCv_tra_mix
	double rhoCv_tra_mix = 0.0;
//#pragma omp parallel for reduction(+:rhoCv_tra_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCv_tra_mix	+= data_Q(s) * species_set.species[s].Cv_tra;
		}
	}


	// 4. rhoCv_rot_mix
	double rhoCv_rot_mix = 0.0;
//#pragma omp parallel for reduction(+:rhoCv_rot_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCv_rot_mix	+= data_Q(s) * species_set.species[s].Cv_rot;
		}
	}


	// 5. rhoCve
	double rhoCve = 0.0;
//#pragma omp parallel for reduction(+:rhoCve)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoCve	= data_Q(s) * species_set.species[s].Cv_tra;
			break;
		}
	}

	double E_h0e = 0.0;
//#pragma omp parallel for reduction(+:E_h0e)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			E_h0e	= data_Q(s) * species_set.species[s].h0;
			break;
		}
	}



	T	= (E - E_h0 - E_k - Er - Eve - Ee) / rhoCv_tra_mix;
	Tr	= Er / rhoCv_rot_mix;
	Tve	= CFD_calculate_Tve(data_Q, species_set, Eve, 10000, 0.0001, CFD_NT);
	Te	= (Ee - E_h0e) / rhoCve;

	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type8: Q => V)");
	}

	if (Te < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Electron Temperature (Type8: Q => V)");
	}

	//Common::ErrorCheckNonNegative<double>(T,   "Temperature(Type8): Q => V");
	//Common::ErrorCheckNonNegative<double>(Tr,  "Rotation Temperature(Type8): Q => V");
	//Common::ErrorCheckNonNegative<double>(Tve, "Vibrational Temperature(Type8): Q => V");
	//Common::ErrorCheckNonNegative<double>(Te,  "Electron Temperature(Type8): Q => V");

	data_V(indexT)		= T;
	data_V(indexTr)		= Tr;
	data_V(indexTve)	= Tve;
	data_V(indexTe)	= Te;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * V ><=> W
 */
void VariableChangeType8::V_to_W(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTr		= indexT + 1;
	int indexTve	= indexT + 2;
	int indexTe		= indexT + 3;

	double T	= data_V(indexT);
	double Tr	= data_V(indexTr);
	double Tve	= data_V(indexTve);
	double Te	= data_V(indexTe);
	double p, pe;


	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoRmix	+= data_V(s)*species_set.species[s].R;
		}
	}

	double rhoeRe	= 0.0;
//#pragma omp parallel for reduction(+:rhoeRe)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoeRe	= data_V(s)*species_set.species[s].R;
			break;
		}
	}


	pe	= rhoeRe * Te;
	p	= rhoRmix*T + pe;

	data_W(indexT)		= p;
	data_W(indexTr)		= Tr;
	data_W(indexTve)	= Tve;
	data_W(indexTe)		= pe;
}

void VariableChangeType8::W_to_V(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTr		= indexT + 1;
	int indexTve	= indexT + 2;
	int indexTe		= indexT + 3;

	double p	= data_W(indexT);
	double Tr	= data_W(indexTr);
	double Tve	= data_W(indexTve);
	double pe	= data_W(indexTe);
	double T, Te;



	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoRmix	+= data_W(s)*species_set.species[s].R;
		}
	}


	double rhoeRe	= 0.0;
//#pragma omp parallel for reduction(+:rhoeRe)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoeRe	= data_W(s)*species_set.species[s].R;
			break;
		}
	}


	Te	= pe / rhoeRe;
	T	= (p - pe) / rhoRmix;

	data_V(indexT)		= T;
	data_V(indexTr)		= Tr;
	data_V(indexTve)	= Tve;
	data_V(indexTe)		= Te;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Q <=> W
 */
void VariableChangeType8::Q_to_W(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTr	= indexT + 1;
	int indexTve	= indexT + 2;
	int indexTe		= indexT + 3;

	double E	= data_Q(indexT);
	double Er	= data_Q(indexTr);
	double Eve	= data_Q(indexTve);
	double Ee	= data_Q(indexTe);

	double T;
	double Tr;
	double Tve;
	double Te;
	double p;
	double pe;



	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_h0	+= data_Q(s)*species_set.species[s].h0;
		}
	}


	// 2. E_k (Kinetic energy)
	double rho_mix	= 0.0;
//#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		E_k	+= data_Q(k)*data_Q(k) / rho_mix;
	}
	E_k *= 0.5;


	// 3. rhoCv_tra_mix
	double rhoCv_tra_mix = 0.0;
//#pragma omp parallel for reduction(+:rhoCv_tra_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCv_tra_mix	+= data_Q(s) * species_set.species[s].Cv_tra;
		}
	}


	// 4. rhoCv_rot_mix
	double rhoCv_rot_mix = 0.0;
//#pragma omp parallel for reduction(+:rhoCv_rot_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCv_rot_mix	+= data_Q(s) * species_set.species[s].Cv_rot;
		}
	}


	// 5. rhoCve
	double rhoCve = 0.0;
//#pragma omp parallel for reduction(+:rhoCve)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoCve	= data_Q(s) * species_set.species[s].Cv_tra;
			break;
		}
	}

	double E_h0e = 0.0;
//#pragma omp parallel for reduction(+:E_h0e)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			E_h0e	= data_Q(s) * species_set.species[s].h0;
			break;
		}
	}



	T	= (E - E_h0 - E_k - Er - Eve - Ee) / rhoCv_tra_mix;
	Tr	= Er / rhoCv_rot_mix;
	Tve	= CFD_calculate_Tve(data_Q, species_set, Eve, 10000, 0.0001, CFD_NT);
	Te	= (Ee - E_h0e) / rhoCve;

	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type8: Q => W)");
	}

	if (Te < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Electron Temperature (Type8: Q => W)");
	}

	//Common::ErrorCheckNonNegative<double>(T,   "Temperature(Type8): Q => V");
	//Common::ErrorCheckNonNegative<double>(Tr,  "Rotation Temperature(Type8): Q => V");
	//Common::ErrorCheckNonNegative<double>(Tve, "Vibrational Temperature(Type8): Q => V");
	//Common::ErrorCheckNonNegative<double>(Te,  "Electron Temperature(Type8): Q => V");


	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoRmix	+= data_Q(s)*species_set.species[s].R;
		}
	}

	double rhoeRe	= 0.0;
//#pragma omp parallel for reduction(+:rhoeRe)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoeRe	= data_Q(s)*species_set.species[s].R;
			break;
		}
	}


	pe	= rhoeRe * Te;
	p	= rhoRmix*T + pe;

	data_W(indexT)		= p;
	data_W(indexTr)		= Tr;
	data_W(indexTve)	= Tve;
	data_W(indexTe)		= pe;
}

void VariableChangeType8::W_to_Q(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTr		= indexT + 1;
	int indexTve	= indexT + 2;
	int indexTe		= indexT + 3;

	double p	= data_W(indexT);
	double Tr	= data_W(indexTr);
	double Tve	= data_W(indexTve);
	double pe	= data_W(indexTe);
	double T, Te;


	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoRmix	+= data_W(s)*species_set.species[s].R;
		}
	}


	double rhoeRe	= 0.0;
//#pragma omp parallel for reduction(+:rhoeRe)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoeRe	= data_W(s)*species_set.species[s].R;
			break;
		}
	}


	Te	= pe / rhoeRe;
	T	= (p - pe) / rhoRmix;


	double E, Er, Eve, Ee;

	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_h0	+= data_W(s)*species_set.species[s].h0;
		}
	}


	// 2. E_k (Kinetic energy)
	double rho_mix	= 0.0;
//#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_W(s);

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		E_k	+= rho_mix*data_W(k)*data_W(k);
	}
	E_k *= 0.5;


	// 3. E_tra (translation/rotation energy) [Except electrons]
	double E_tra = 0.0;
//#pragma omp parallel for reduction(+:E_tra)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_tra	+= data_W(s) * (species_set.species[s].Cv_tra*T);
		}
	}


	// 4. E_rot (translation/rotation energy) [Except electrons]
	Er = 0.0;
//#pragma omp parallel for reduction(+:Er)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			Er	+= data_W(s) * (species_set.species[s].Cv_rot*Tr);
		}
	}


	// 5. E_VE (VIBRATION-ELECTRONIC energy)
	Eve = 0.0;
//#pragma omp parallel for reduction(+:Eve)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			Eve	+= data_W(s) * species_set.species[s].e_VE(Tve);
		}
	}


	// 5. Ee
	Ee = 0.0;
//#pragma omp parallel for reduction(+:Ee)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			Ee	= data_W(s) * (species_set.species[s].Cv_tra*Te + species_set.species[s].h0);
			break;
		}
	}

	E	= E_tra + E_h0 + E_k + Er + Eve + Ee;

	data_Q(indexT)		= E;
	data_Q(indexTr)		= Er;
	data_Q(indexTve)	= Eve;
	data_Q(indexTe)		= Ee;
}



void VariableChangeType8::From_Q(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
{
	// 1. Calculate rho
	double rho = 0.0;
//#pragma omp parallel for reduction(+:rho) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rho	+= data_Q(s);


	// 2. Calculate rhoCv and rhoR
	double rhoCvtra_wo_e	= 0.0;
	double rhoCvrot_wo_e	= 0.0;
	double rhoCve 			= 0.0;
	double rhoR_wo_e		= 0.0;
	double rhoRe			= 0.0;
	double E_h0				= 0.0;
	double E_h0e			= 0.0;

//#pragma omp parallel for reduction(+:rhoCvtra_wo_e, rhoCvrot_wo_e, rhoCve, rhoR_wo_e, rhoRe) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCvtra_wo_e	+= data_Q(s)*species_set.species[s].Cv_tra;
			rhoCvrot_wo_e	+= data_Q(s)*species_set.species[s].Cv_rot;
			rhoR_wo_e	+= data_Q(s)*species_set.species[s].R;

			E_h0		+= data_Q(s)*species_set.species[s].h0;
		}
		else
		{
			rhoCve	+= data_Q(s)*species_set.species[s].Cv_tra;
			rhoRe	+= data_Q(s)*species_set.species[s].R;

			E_h0e	+= data_Q(s)*species_set.species[s].h0;
		}
	}



	// 1. Species and velocities
//#pragma omp parallel for num_threads(CFD_NT)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_V(s)	= data_Q(s);
		data_W(s)	= data_Q(s);
	}

//#pragma omp parallel for num_threads(CFD_NT)
#pragma ivdep
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		data_V(k)	= data_Q(k) / rho;
		data_W(k)	= data_V(k);
	}

	// 2. Calculate E_ho and E_k
	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_V(k)*data_V(k);
	E_k = 0.5*rho*E_k;

	// 3. Calculate T and Te
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;
	int indexTv	= species_set.NS + ND + 2;
	int indexTe	= species_set.NS + ND + 3;

	double T	= (data_Q(indexT) - data_Q(indexTr) - data_Q(indexTv) - data_Q(indexTe) - E_k - E_h0) / rhoCvtra_wo_e;
	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type8: fromQ => ALL)");
	}
	//Common::ErrorCheckNonNegative<double>(T, "Temperature(Type8)");

	double Tr	= data_Q(indexTr) / rhoCvrot_wo_e;
	double Tv	= CFD_calculate_Tve(data_Q, species_set, data_Q(indexTv), 10000, 0.0001, CFD_NT);
	double Te	= (data_Q(indexTe) - E_h0e) / rhoCve;

	data_V(indexT)	= T;
	data_V(indexTr)	= Tr;
	data_V(indexTv)	= Tv;
	data_V(indexTe)	= Te;

	double pe		= rhoRe*Te;
	data_W(indexT)	= rhoR_wo_e*T + pe;
	data_W(indexTr)	= Tr;
	data_W(indexTv)	= Tv;
	data_W(indexTe)	= pe;


	// Assign Mixture properties
	// Calculate N
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Xs(s)	=	data_Q(s) / species_set.species[s].m;
	}

	double N = 0.0;
#pragma omp parallel for reduction(+:N) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		N	+= data_Xs(s);
	}

//#pragma omp parallel for num_threads(CFD_NT)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Xs(s)	/= N;
		data_Ys(s)	= data_Q(s) / rho;
	}


	data_MIX(0)	= rho;
	data_MIX(1) = rhoR_wo_e;
	data_MIX(2) = rho / N * AMU_SI;
	data_MIX(3) = 1.0 + rhoR_wo_e/(rhoCvtra_wo_e + rhoCvrot_wo_e);
	data_MIX(4) = rhoCvtra_wo_e;
	data_MIX(5) = rhoCvrot_wo_e;
}



void VariableChangeType8::From_V(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
{
	// 1. Calculate rho
	double rho = 0.0;
//#pragma omp parallel for reduction(+:rho) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rho	+= data_V(s);


	// 2. Calculate rhoCv and rhoR
	double rhoCvtra_wo_e	= 0.0;
	double rhoCvrot_wo_e	= 0.0;
	double rhoCve 			= 0.0;
	double rhoR_wo_e		= 0.0;
	double rhoRe			= 0.0;
	double E_h0				= 0.0;
	double E_h0e			= 0.0;

//#pragma omp parallel for reduction(+:rhoCvtra_wo_e, rhoCvrot_wo_e, rhoCve, rhoR_wo_e, rhoRe) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCvtra_wo_e	+= data_V(s)*species_set.species[s].Cv_tra;
			rhoCvrot_wo_e	+= data_V(s)*species_set.species[s].Cv_rot;
			rhoR_wo_e		+= data_V(s)*species_set.species[s].R;

			E_h0			+= data_V(s)*species_set.species[s].h0;
		}
		else
		{
			rhoCve	+= data_V(s)*species_set.species[s].Cv_tra;
			rhoRe	+= data_V(s)*species_set.species[s].R;

			E_h0e	+= data_V(s)*species_set.species[s].h0;
		}
	}



	// 1. Species and velocities
//#pragma omp parallel for num_threads(CFD_NT)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Q(s)	= data_V(s);
		data_W(s)	= data_V(s);
	}

//#pragma omp parallel for num_threads(CFD_NT)
#pragma ivdep
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		data_Q(k)	= data_V(k) * rho;
		data_W(k)	= data_V(k);
	}

	// 2. Calculate E_ho and E_k
	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_V(k)*data_V(k);
	E_k = 0.5*rho*E_k;

	// 3. Calculate T and Te
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;
	int indexTv	= species_set.NS + ND + 2;
	int indexTe	= species_set.NS + ND + 3;

	double T	= data_V(indexT);
	double Tr	= data_V(indexTr);
	double Tv	= data_V(indexTv);
	double Te	= data_V(indexTe);


	double Er = rhoCvrot_wo_e	* Tr;
	double Ev = 0.0;
//#pragma omp parallel for reduction(+:Ev) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			Ev	+= data_V(s) * species_set.species[s].e_VE(Tv);
		}
	}
	double Ee = rhoCve*Te + E_h0e;
	data_Q(indexT)	= rhoCvtra_wo_e*T + E_h0 + E_k + Er + Ev + Ee;
	data_Q(indexTr)	= Er;
	data_Q(indexTv)	= Ev;
	data_Q(indexTe)	= Ee;

	double pe		= rhoRe*Te;
	data_W(indexT)	= rhoR_wo_e*T + pe;
	data_W(indexTr)	= Tr;
	data_W(indexTv)	= Tv;
	data_W(indexTe)	= pe;


	// Assign Mixture properties
	// Calculate N
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Xs(s)	=	data_Q(s) / species_set.species[s].m;
	}

	double N = 0.0;
#pragma omp parallel for reduction(+:N) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		N	+= data_Xs(s);
	}

//#pragma omp parallel for num_threads(CFD_NT)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Xs(s)	/= N;
		data_Ys(s)	= data_Q(s) / rho;
	}


	data_MIX(0)	= rho;
	data_MIX(1) = rhoR_wo_e;
	data_MIX(2) = rho / N * AMU_SI;
	data_MIX(3) = 1.0 + rhoR_wo_e/(rhoCvtra_wo_e + rhoCvrot_wo_e);
	data_MIX(4) = rhoCvtra_wo_e;
	data_MIX(5) = rhoCvrot_wo_e;
}


void VariableChangeType8::From_W(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
{
	// 1. Calculate rho
	double rho = 0.0;
//#pragma omp parallel for reduction(+:rho) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rho	+= data_W(s);


	// 2. Calculate rhoCv and rhoR
	double rhoCvtra_wo_e	= 0.0;
	double rhoCvrot_wo_e	= 0.0;
	double rhoCve 			= 0.0;
	double rhoR_wo_e		= 0.0;
	double rhoRe			= 0.0;
	double E_h0				= 0.0;
	double E_h0e			= 0.0;

//#pragma omp parallel for reduction(+:rhoCvtra_wo_e, rhoCvrot_wo_e, rhoCve, rhoR_wo_e, rhoRe) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCvtra_wo_e	+= data_W(s)*species_set.species[s].Cv_tra;
			rhoCvrot_wo_e	+= data_W(s)*species_set.species[s].Cv_rot;
			rhoR_wo_e		+= data_W(s)*species_set.species[s].R;

			E_h0			+= data_W(s)*species_set.species[s].h0;
		}
		else
		{
			rhoCve	+= data_W(s)*species_set.species[s].Cv_tra;
			rhoRe	+= data_W(s)*species_set.species[s].R;

			E_h0e	+= data_W(s)*species_set.species[s].h0;
		}
	}



	// 1. Species and velocities
//#pragma omp parallel for num_threads(CFD_NT)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Q(s)	= data_W(s);
		data_V(s)	= data_W(s);
	}

//#pragma omp parallel for num_threads(CFD_NT)
#pragma ivdep
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		data_Q(k)	= data_W(k) * rho;
		data_V(k)	= data_W(k);
	}

	// 2. Calculate E_ho and E_k
	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_V(k)*data_V(k);
	E_k = 0.5*rho*E_k;

	// 3. Calculate T and Te
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;
	int indexTv	= species_set.NS + ND + 2;
	int indexTe	= species_set.NS + ND + 3;

	double T	= (data_W(indexT) - data_W(indexTe)) / rhoR_wo_e;
	double Tr	= data_W(indexTr);
	double Tv	= data_W(indexTv);
	double Te	= data_W(indexTe) / rhoRe;


	double Er = rhoCvrot_wo_e	* Tr;
	double Ev = 0.0;
//#pragma omp parallel for reduction(+:Ev) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			Ev	+= data_V(s) * species_set.species[s].e_VE(Tv);
		}
	}
	double Ee = rhoCve*Te + E_h0e;
	data_Q(indexT)	= rhoCvtra_wo_e*T + E_h0 + E_k + Er + Ev + Ee;
	data_Q(indexTr)	= Er;
	data_Q(indexTv)	= Ev;
	data_Q(indexTe)	= Ee;

	data_V(indexT)	= T;
	data_V(indexTr)	= Tr;
	data_V(indexTv)	= Tv;
	data_V(indexTe)	= Te;


	// Assign Mixture properties
	// Calculate N
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Xs(s)	=	data_Q(s) / species_set.species[s].m;
	}

	double N = 0.0;
#pragma omp parallel for reduction(+:N) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		N	+= data_Xs(s);
	}

//#pragma omp parallel for num_threads(CFD_NT)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Xs(s)	/= N;
		data_Ys(s)	= data_Q(s) / rho;
	}


	data_MIX(0)	= rho;
	data_MIX(1) = rhoR_wo_e;
	data_MIX(2) = rho / N * AMU_SI;
	data_MIX(3) = 1.0 + rhoR_wo_e/(rhoCvtra_wo_e + rhoCvrot_wo_e);
	data_MIX(4) = rhoCvtra_wo_e;
	data_MIX(5) = rhoCvrot_wo_e;
}

}
}
