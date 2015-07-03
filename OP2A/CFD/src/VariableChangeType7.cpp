/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 19, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChangeType7.cpp
 * 			-  
 *  
 */


#include <omp.h>

#include "CFD/include/VariableChangeTypes.hpp"


namespace OP2A{
namespace CFD{


void VariableChangeType7::V_to_Q(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTr		= species_set.NS + ND + 1;
	int indexTve	= species_set.NS + ND + 2;

	double T	= data_V(indexT);
	double Tr	= data_V(indexTr);
	double Tve	= data_V(indexTve);

	double E, Er, Eve;



	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_h0	+= data_V(s)*species_set.species[s].h0;
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


	// 3. E_tra (translation/rotation energy)
	double E_tra = 0.0;
//#pragma omp parallel for reduction(+:E_tra)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_tra	+= data_V(s) * (species_set.species[s].Cv_tra*T);
	}


	// 4. Er (Rotation energy)
	Er = 0.0;
//#pragma omp parallel for reduction(+:Er)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Er	+= data_V(s) * (species_set.species[s].Cv_rot*Tr);
	}


	// 5. E_VE (VIBRATION-ELECTRONIC energy)
	Eve = 0.0;
//#pragma omp parallel for reduction(+:Eve)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Eve	+= data_V(s) * species_set.species[s].e_VE(Tve);
	}

	E	= E_tra + E_k + E_h0 + Er + Eve;
	data_Q(indexT)		= E;
	data_Q(indexTr)		= Er;
	data_Q(indexTve)	= Eve;
}

void VariableChangeType7::Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTr		= species_set.NS + ND + 1;
	int indexTve	= species_set.NS + ND + 2;

	double E	= data_Q(indexT);
	double Er	= data_Q(indexTr);
	double Eve	= data_Q(indexTve);
	double T, Tr, Tve;


	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_h0	+= data_Q(s)*species_set.species[s].h0;
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
		rhoCv_tra_mix	+= data_Q(s) * species_set.species[s].Cv_tra;
	}


	// 4. rhoCv_rot_mix
	double rhoCv_rot_mix = 0.0;
//#pragma omp parallel for reduction(+:rhoCv_rot_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		rhoCv_rot_mix	+= data_Q(s) * species_set.species[s].Cv_rot;
	}


	T	= (E - Er - Eve - E_h0 - E_k) / rhoCv_tra_mix;
	Tr	= Er / rhoCv_rot_mix;
	Tve	= CFD_calculate_Tve(data_Q, species_set, Eve, 10000, 0.01, CFD_NT);

	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type7: Q => V)");
	}
	//Common::ErrorCheckNonNegative<double>(T,   "Temperature(Type7): Q => V");
	//Common::ErrorCheckNonNegative<double>(Tr,  "Rotational Temperature(Type7): Q => V");
	//Common::ErrorCheckNonNegative<double>(Tve, "Vibrational Temperature(Type7): Q => V");

	data_V(indexT)		= T;
	data_V(indexTr)		= Tr;
	data_V(indexTve)	= Tve;
}



//////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * V ><=> W
 */
void VariableChangeType7::V_to_W(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTr		= species_set.NS + ND + 1;
	int indexTve	= species_set.NS + ND + 2;

	double T	= data_V(indexT);
	double Tr	= data_V(indexTr);
	double Tve	= data_V(indexTve);

	double p;


	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_V(s)*species_set.species[s].R;

	p	= rhoRmix*T;

	data_W(indexT)		= p;
	data_W(indexTr)		= Tr;
	data_W(indexTve)	= Tve;
}

void VariableChangeType7::W_to_V(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTr		= species_set.NS + ND + 1;
	int indexTve	= species_set.NS + ND + 2;

	double p	= data_W(indexT);
	double Tr	= data_W(indexTr);
	double Tve	= data_W(indexTve);
	double T;


	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_V(s)*species_set.species[s].R;

	T	= p / rhoRmix;

	data_V(indexT)		= T;
	data_V(indexTr)		= Tr;
	data_V(indexTve)	= Tve;
}





/////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Q <=> W
 */
void VariableChangeType7::Q_to_W(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTr		= species_set.NS + ND + 1;
	int indexTve	= species_set.NS + ND + 2;

	double E	= data_Q(indexT);
	double Er	= data_Q(indexTr);
	double Eve	= data_Q(indexTve);
	double T, Tr, Tve;


	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_h0	+= data_Q(s)*species_set.species[s].h0;
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
		rhoCv_tra_mix	+= data_Q(s) * species_set.species[s].Cv_tra;
	}


	// 4. rhoCv_rot_mix
	double rhoCv_rot_mix = 0.0;
//#pragma omp parallel for reduction(+:rhoCv_rot_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		rhoCv_rot_mix	+= data_Q(s) * species_set.species[s].Cv_rot;
	}


	T	= (E - Er - Eve - E_h0 - E_k) / rhoCv_tra_mix;
	Tr	= Er / rhoCv_rot_mix;
	Tve	= CFD_calculate_Tve(data_Q, species_set, Eve, 10000, 0.0001, CFD_NT);

	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type7: Q => W)");
	}

	//Common::ErrorCheckNonNegative<double>(T,   "Temperature(Type7): Q => V");
	//Common::ErrorCheckNonNegative<double>(Tr,  "Rotational Temperature(Type7): Q => V");
	//Common::ErrorCheckNonNegative<double>(Tve, "Vibrational Temperature(Type7): Q => V");


	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_Q(s)*species_set.species[s].R;

	double p	= rhoRmix*T;

	data_W(indexT)		= p;
	data_W(indexTr)		= Tr;
	data_W(indexTve)	= Tve;
}



void VariableChangeType7::W_to_Q(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTr		= species_set.NS + ND + 1;
	int indexTve	= species_set.NS + ND + 2;

	double p	= data_W(indexT);
	double Tr	= data_W(indexTr);
	double Tve	= data_W(indexTve);
	double T;
	double E, Er, Eve;


	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_W(s)*species_set.species[s].R;

	T	= p / rhoRmix;


	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_h0	+= data_W(s)*species_set.species[s].h0;
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


	// 3. E_tra (translation/rotation energy)
	double E_tra = 0.0;
//#pragma omp parallel for reduction(+:E_tra)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_tra	+= data_W(s) * (species_set.species[s].Cv_tra*T);
	}


	// 4. Er (Rotation energy)
	Er = 0.0;
//#pragma omp parallel for reduction(+:Er)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Er	+= data_W(s) * (species_set.species[s].Cv_rot*Tr);
	}


	// 5. E_VE (VIBRATION-ELECTRONIC energy)
	Eve = 0.0;
//#pragma omp parallel for reduction(+:Eve)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Eve	+= data_W(s) * species_set.species[s].e_VE(Tve);
	}

	E	= E_tra + E_k + E_h0 + Er + Eve;
	data_Q(indexT)		= E;
	data_Q(indexTr)		= Er;
	data_Q(indexTve)	= Eve;
}



void VariableChangeType7::From_Q(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
{
	// 1. Calculate rho
	double rho = 0.0;
//#pragma omp parallel for reduction(+:rho) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rho	+= data_Q(s);


	// 2. Calculate rhoCv and rhoR
	double rhoCvtra = 0.0;
	double rhoCvrot = 0.0;
	double rhoR  = 0.0;

//#pragma omp parallel for reduction(+:rhoCvtra) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoCvtra	+= data_Q(s)*species_set.species[s].Cv_tra;

//#pragma omp parallel for reduction(+:rhoCvrot) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoCvrot	+= data_Q(s)*species_set.species[s].Cv_rot;

//#pragma omp parallel for reduction(+:rhoR) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoR	+= data_Q(s)*species_set.species[s].R;



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
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_Q(s)*species_set.species[s].h0;

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_V(k)*data_V(k);
	E_k = 0.5*rho*E_k;


	// 3. Calculate T
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;
	int indexTv	= species_set.NS + ND + 2;

	double T	= (data_Q(indexT) - data_Q(indexTr) - data_Q(indexTv) - E_k - E_h0) / rhoCvtra;
	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type7: fromQ => ALL)");
	}
	//Common::ErrorCheckNonNegative<double>(T, "Temperature(Type1)");

	double Tr	= data_Q(indexTr) / rhoCvrot;
	double Tv	= CFD_calculate_Tve(data_Q, species_set, data_Q(indexTv), 10000, 0.01, CFD_NT);


	data_V(indexT)	= T;
	data_V(indexTr)	= Tr;
	data_V(indexTv)	= Tv;

	data_W(indexT)	= rhoR*T;
	data_W(indexTr)	= Tr;
	data_W(indexTv)	= Tv;



	// Assign Mixture properties
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Xs(s)	=	data_Q(s) / species_set.species[s].m;
	}

	// Calculate N
	double N = 0.0;
//#pragma omp parallel for reduction(+:N) num_threads(CFD_NT)
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
	data_MIX(1) = rhoR;
	data_MIX(2) = rho / N * AMU_SI;
	data_MIX(3) = 1.0 + rhoR/(rhoCvtra + rhoCvrot);
	data_MIX(4) = rhoCvtra;
	data_MIX(5) = rhoCvrot;
}


void VariableChangeType7::From_V(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
{
	// 1. Calculate rho
	double rho = 0.0;
//#pragma omp parallel for reduction(+:rho) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rho	+= data_V(s);


	// 2. Calculate rhoCv and rhoR
	double rhoCvtra = 0.0;
	double rhoCvrot = 0.0;
	double rhoR  = 0.0;

//#pragma omp parallel for reduction(+:rhoCvtra) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoCvtra	+= data_V(s)*species_set.species[s].Cv_tra;

//#pragma omp parallel for reduction(+:rhoCvrot) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoCvrot	+= data_V(s)*species_set.species[s].Cv_rot;

//#pragma omp parallel for reduction(+:rhoR) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoR	+= data_V(s)*species_set.species[s].R;



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
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_Q(s)*species_set.species[s].h0;

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_V(k)*data_V(k);
	E_k = 0.5*rho*E_k;


	// 3. Calculate T
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;
	int indexTv	= species_set.NS + ND + 2;

	double T	= data_V(indexT);
	double Tr	= data_V(indexTr);
	double Tv	= data_V(indexTv);

	double Er	= rhoCvrot * Tr;
	double Ev	= 0.0;
//#pragma omp parallel for reduction(+:Ev) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Ev	+= data_V(s) * species_set.species[s].e_VE(Tv);
	}

	data_Q(indexT)	= rhoCvtra*T + Er + Ev + E_h0 + E_k;
	data_Q(indexTr)	= Er;
	data_Q(indexTv)	= Ev;

	data_W(indexT)	= rhoR*T;
	data_W(indexTr)	= Tr;
	data_W(indexTv)	= Tv;



	// Assign Mixture properties
	// Calculate N
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Xs(s)	=	data_Q(s) / species_set.species[s].m;
	}


	double N = 0.0;
//#pragma omp parallel for reduction(+:N) num_threads(CFD_NT)
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
	data_MIX(1) = rhoR;
	data_MIX(2) = rho / N * AMU_SI;
	data_MIX(3) = 1.0 + rhoR/(rhoCvtra + rhoCvrot);
	data_MIX(4) = rhoCvtra;
	data_MIX(5) = rhoCvrot;
}


void VariableChangeType7::From_W(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
{
	// 1. Calculate rho
	double rho = 0.0;
//#pragma omp parallel for reduction(+:rho) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rho	+= data_W(s);


	// 2. Calculate rhoCv and rhoR
	double rhoCvtra = 0.0;
	double rhoCvrot = 0.0;
	double rhoR  = 0.0;

//#pragma omp parallel for reduction(+:rhoCvtra) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoCvtra	+= data_W(s)*species_set.species[s].Cv_tra;

//#pragma omp parallel for reduction(+:rhoCvrot) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoCvrot	+= data_W(s)*species_set.species[s].Cv_rot;

//#pragma omp parallel for reduction(+:rhoR) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoR	+= data_W(s)*species_set.species[s].R;



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
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_Q(s)*species_set.species[s].h0;

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_V(k)*data_V(k);
	E_k = 0.5*rho*E_k;


	// 3. Calculate T
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;
	int indexTv	= species_set.NS + ND + 2;

	double T	= data_W(indexT) / rhoR;
	double Tr	= data_W(indexTr);
	double Tv	= data_W(indexTv);

	double Er	= rhoCvrot * Tr;
	double Ev	= 0.0;
//#pragma omp parallel for reduction(+:Ev) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Ev	+= data_W(s) * species_set.species[s].e_VE(Tv);
	}

	data_Q(indexT)	= rhoCvtra*T + Er + Ev + E_h0 + E_k;
	data_Q(indexTr)	= Er;
	data_Q(indexTv)	= Ev;

	data_V(indexT)	= T;
	data_V(indexTr)	= Tr;
	data_V(indexTv)	= Tv;



	// Assign Mixture properties
	// Calculate N
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Xs(s)	=	data_Q(s) / species_set.species[s].m;
	}

	double N = 0.0;
//#pragma omp parallel for reduction(+:N) num_threads(CFD_NT)
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
	data_MIX(1) = rhoR;
	data_MIX(2) = rho / N * AMU_SI;
	data_MIX(3) = 1.0 + rhoR/(rhoCvtra + rhoCvrot);
	data_MIX(4) = rhoCvtra;
	data_MIX(5) = rhoCvrot;
}



}
}

