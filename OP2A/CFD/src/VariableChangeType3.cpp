/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 19, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChangeType3.cpp
 * 			-  
 *  
 */



#include <omp.h>

#include "CFD/include/VariableChangeTypes.hpp"


namespace OP2A{
namespace CFD{



void VariableChangeType3::V_to_Q(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	int indexT		= species_set.NS + ND;
	int indexTve	= indexT + 1;

	double T	= data_V(indexT);
	double Tve	= data_V(indexTve);


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


	// 3. E_int (translation/rotation energy) [Except electrons]
	double E_int = 0.0;
//#pragma omp parallel for reduction(+:E_int)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_int	+= data_V(s) * (species_set.species[s].Cv_tr*T);
	}


	// 4. E_VE (VIBRATION-ELECTRONIC energy)
	double E_ve = 0.0;
//#pragma omp parallel for reduction(+:E_ve)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_ve	+= data_V(s) * species_set.species[s].e_VE(Tve);
	}


	data_Q(indexT)	= E_int + E_k + E_h0 + E_ve;
	data_Q(indexTve)	= E_ve;
}

void VariableChangeType3::Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	double rho_mix	= 0.0;
//#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);


	int indexT		= species_set.NS + ND;
	int indexTve	= indexT + 1;

	double E	= data_Q(indexT);
	double Eve  = data_Q(indexTve);

	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_h0	+= data_Q(s)*species_set.species[s].h0;
	}


	// 2. E_k (Kinetic energy)
	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		E_k	+= data_Q(k)*data_Q(k) / rho_mix;
	}
	E_k *= 0.5;


	// 3. rhoCvmix
	double rhoCvmix = 0.0;
//#pragma omp parallel for reduction(+:rhoCvmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		rhoCvmix	+= data_Q(s) * species_set.species[s].Cv_tr;
	}


	double T, Tve;
	T	= (E - E_h0 - E_k - Eve) / rhoCvmix;
	Tve	= CFD_calculate_Tve(data_Q, species_set, Eve, 10000, 0.0001, CFD_NT);

	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type3: Q => V)");
	}

	//Common::ErrorCheckNonNegative<double>(T,   "Temperature(Type3): Q => V");
	//Common::ErrorCheckNonNegative<double>(Tve, "Vibrational Temperature(Type3): Q => V");

	data_V(indexT)		= T;
	data_V(indexTve)	= Tve;
}



//////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * V ><=> W
 */
void VariableChangeType3::V_to_W(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_V(s)*species_set.species[s].R;

	int indexT		= species_set.NS + ND;
	int indexTve	= species_set.NS + ND + 1;

	double T	= data_V(indexT);
	double Tve	= data_V(indexTve);
	double p	= rhoRmix*T;

	data_W(indexT)	= p;
	data_W(indexTve)	= Tve;
}

void VariableChangeType3::W_to_V(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_V(s)*species_set.species[s].R;


	int indexT		= species_set.NS + ND;
	int indexTve	= species_set.NS + ND + 1;

	double p	= data_W(indexT);
	double Tve	= data_W(indexTve);
	double T	= p / rhoRmix;

	data_V(indexT)		= T;
	data_V(indexTve)	= Tve;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Q <=> W
 */
void VariableChangeType3::Q_to_W(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{


	int indexT		= species_set.NS + ND;
	int indexTve	= indexT + 1;

	double E	= data_Q(indexT);
	double Eve  = data_Q(indexTve);



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


	// 3. rhoCvmix
	double rhoCvmix = 0.0;
//#pragma omp parallel for reduction(+:rhoCvmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		rhoCvmix	+= data_Q(s) * species_set.species[s].Cv_tr;
	}



	double p;
	double T;
	double Tve;

	T	= (E - E_h0 - E_k - Eve) / rhoCvmix;
	Tve	= CFD_calculate_Tve(data_Q, species_set, Eve, 10000, 0.0001, CFD_NT);

	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type3: Q => W)");
	}

	//Common::ErrorCheckNonNegative<double>(T,   "Temperature(Type3): Q => W");
	//Common::ErrorCheckNonNegative<double>(Tve, "Vibrational Temperature(Type3): Q => Q");


	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_Q(s)*species_set.species[s].R;

	p	= rhoRmix*T;

	data_W(indexT)		= p;
	data_W(indexTve)	= Tve;
}


void VariableChangeType3::W_to_Q(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_W(s)*species_set.species[s].R;

	int indexT		= species_set.NS + ND;
	int indexTve	= indexT + 1;

	double p	= data_W(indexT);
	double Tve	= data_W(indexTve);
	double T	= p / rhoRmix;


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


	// 3. E_int (translation/rotation energy) [Except electrons]
	double E_int = 0.0;
//#pragma omp parallel for reduction(+:E_int)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_int	+= data_W(s) * (species_set.species[s].Cv_tr*T);
	}


	// 4. E_VE (VIBRATION-ELECTRONIC energy)
	double E_ve = 0.0;
//#pragma omp parallel for reduction(+:E_ve)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_ve	+= data_W(s) * species_set.species[s].e_VE(Tve);
	}


	data_Q(indexT)	= E_int + E_k + E_h0 + E_ve;
	data_Q(indexT)	= E_ve;
}


void VariableChangeType3::From_Q(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
{
	// 1. Calculate rho
	double rho = 0.0;
//#pragma omp parallel for reduction(+:rho) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rho	+= data_Q(s);


	// 2. Calculate rhoCv and rhoR
	double rhoCv = 0.0;
	double rhoR  = 0.0;

//#pragma omp parallel for reduction(+:rhoCv) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoCv	+= data_Q(s)*species_set.species[s].Cv_tr;

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
	int indexT		= species_set.NS + ND;
	int indexTve	= species_set.NS + ND + 1;

	double T	= (data_Q(indexT) - data_Q(indexTve) - E_k - E_h0) / rhoCv;
	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type3: fromQ => ALL)");
	}
	//Common::ErrorCheckNonNegative<double>(T, "Temperature(Type1)");

	double Tve	= CFD_calculate_Tve(data_Q, species_set, data_Q(indexTve), 10000, 0.01, CFD_NT);


	data_V(indexT)		= T;
	data_V(indexTve)	= Tve;

	data_W(indexT)		= rhoR*T;
	data_W(indexTve)	= Tve;



	// Assign Mixture properties
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Xs(s) = data_Q(s) / species_set.species[s].m;
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
	data_MIX(3) = 1.0 + rhoR/rhoCv;
	data_MIX(4) = rhoCv;
}


void VariableChangeType3::From_V(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
{
	// 1. Calculate rho
	double rho = 0.0;
//#pragma omp parallel for reduction(+:rho) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rho	+= data_V(s);


	// 2. Calculate rhoCv and rhoR
	double rhoCv = 0.0;
	double rhoR  = 0.0;

//#pragma omp parallel for reduction(+:rhoCv) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoCv	+= data_V(s)*species_set.species[s].Cv_tr;

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


	// 3. Calculate E /Eve
	int indexT		= species_set.NS + ND;
	int indexTve	= species_set.NS + ND + 1;

	double T	= data_V(indexT);
	double Tve	= data_V(indexTve);

	double Eve	= 0.0;
//#pragma omp parallel for reduction(+:Eve) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Eve	+= data_V(s) * species_set.species[s].e_VE(Tve);
	}


	data_Q(indexT)		= rhoCv*T + E_k + E_h0 + Eve;
	data_Q(indexTve)	= Eve;

	data_W(indexT)		= rhoR*T;
	data_W(indexTve)	= Tve;



	// Assign Mixture properties
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{

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
	data_MIX(3) = 1.0 + rhoR/rhoCv;
	data_MIX(4) = rhoCv;
}


void VariableChangeType3::From_W(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
{
	// 1. Calculate rho
	double rho = 0.0;
//#pragma omp parallel for reduction(+:rho) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rho	+= data_W(s);


	// 2. Calculate rhoCv and rhoR
	double rhoCv = 0.0;
	double rhoR  = 0.0;

//#pragma omp parallel for reduction(+:rhoCv) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoCv	+= data_W(s)*species_set.species[s].Cv_tr;

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


	// 3. Calculate E /Eve
	int indexT		= species_set.NS + ND;
	int indexTve	= species_set.NS + ND + 1;

	double T	= data_W(indexT) / rhoR;
	double Tve	= data_W(indexTve);

	double Eve	= 0.0;
//#pragma omp parallel for reduction(+:Eve) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Eve	+= data_W(s) * species_set.species[s].e_VE(Tve);
	}


	data_Q(indexT)		= rhoCv*T + E_k + E_h0 + Eve;
	data_Q(indexTve)	= Eve;

	data_V(indexT)		= T;
	data_V(indexTve)	= Tve;



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
	data_MIX(3) = 1.0 + rhoR/rhoCv;
	data_MIX(4) = rhoCv;
}



}
}
