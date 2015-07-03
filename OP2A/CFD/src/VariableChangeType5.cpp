/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 19, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChangeType5.cpp
 * 			-  
 *  
 */



#include <omp.h>

#include "CFD/include/VariableChangeTypes.hpp"


namespace OP2A{
namespace CFD{




void VariableChangeType5::V_to_Q(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;

	double T	= data_V(indexT);
	double Tr	= data_V(indexTr);

	double E, Er;


	// 1. Enthalpy of formation
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_V(s)*species_set.species[s].h0;


	// 2. Kinetic energy
	double rho_mix	= 0.0;
//#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_V(s);

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= rho_mix*data_V(k)*data_V(k);
	E_k *= 0.5;


	// 3. Translation
	double E_tra = 0.0;
//#pragma omp parallel for reduction(+:E_tra)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_tra	+= data_V(s) * (species_set.species[s].Cv_tra*T);
	}


	// 4. Rotation
	Er = 0.0;
//#pragma omp parallel for reduction(+:Er)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Er	+= data_V(s) * (species_set.species[s].Cv_rot*Tr);
	}


	E = E_tra + E_h0 + E_k + Er;

	data_Q(indexT)	= E;
	data_Q(indexTr)	= Er;
}


void VariableChangeType5::Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;

	double E	= data_Q(indexT);
	double Er	= data_Q(indexTr);

	double T, Tr;



	// Enthalpy of formation
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_Q(s)*species_set.species[s].h0;


	// Kinetic Energy
	double rho_mix	= 0.0;
//#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_Q(k)*data_Q(k) / rho_mix;
	E_k *= 0.5;


	// 1. Cv_tra_mix
	double Cv_tra_mix = 0.0;
//#pragma omp parallel for reduction(+:Cv_tra_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Cv_tra_mix	+= data_Q(s) * species_set.species[s].Cv_tra;
	}

	// 2. Cv_rot_mix
	double Cv_rot_mix = 0.0;
//#pragma omp parallel for reduction(+:Cv_rot_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Cv_rot_mix	+= data_Q(s) * species_set.species[s].Cv_rot;
	}


	T	= (E - E_k - E_h0 - Er) / Cv_tra_mix;
	Tr	= Er / Cv_rot_mix;

	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type5: Q => V)");
	}

	//Common::ErrorCheckNonNegative<double>(T, "Temperature(Type5): Q => V");
	//Common::ErrorCheckNonNegative<double>(Tr, "Rotational Temperature(Type5): Q => V");

	data_V(indexT)	= T;
	data_V(indexTr)	= Tr;
}




//////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * V ><=> W
 */

void VariableChangeType5::V_to_W(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;

	double T	= data_V(indexT);
	double Tr	= data_V(indexTr);

	double p;


	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_V(s)*species_set.species[s].R;

	p	= rhoRmix * T;

	// Pressure
	data_W(indexT)	= p;
	data_W(indexTr)	= Tr;
}

void VariableChangeType5::W_to_V(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;

	double p	= data_W(indexT);
	double Tr	= data_W(indexTr);
	double T;


	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_W(s)*species_set.species[s].R;
	T	= p / rhoRmix;

	data_V(indexT)	= T;
	data_V(indexTr)	= Tr;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Q <=> W
 */
void VariableChangeType5::Q_to_W(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;

	double E	= data_Q(indexT);
	double Er	= data_Q(indexTr);

	double T, Tr, p;



	// Enthalpy of formation
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_Q(s)*species_set.species[s].h0;


	// Kinetic Energy
	double rho_mix	= 0.0;
//#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_Q(k)*data_Q(k) / rho_mix;
	E_k *= 0.5;


	// 1. Cv_tra_mix
	double Cv_tra_mix = 0.0;
//#pragma omp parallel for reduction(+:Cv_tra_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Cv_tra_mix	+= data_Q(s) * species_set.species[s].Cv_tra;
	}

	// 2. Cv_rot_mix
	double Cv_rot_mix = 0.0;
//#pragma omp parallel for reduction(+:Cv_rot_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Cv_rot_mix	+= data_Q(s) * species_set.species[s].Cv_rot;
	}


	T	= (E - E_k - E_h0 - Er) / Cv_tra_mix;
	Tr	= Er / Cv_rot_mix;

	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type5: Q => W)");
	}
	//Common::ErrorCheckNonNegative<double>(T, "Temperature(Type5): Q => V");
	//Common::ErrorCheckNonNegative<double>(Tr, "Rotational Temperature(Type5): Q => V");

	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_Q(s)*species_set.species[s].R;
	p	= rhoRmix * T;

	data_W(indexT)	= p;
	data_W(indexTr)	= Tr;
}



void VariableChangeType5::W_to_Q(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;

	double p	= data_W(indexT);
	double Tr	= data_W(indexTr);
	double T;


	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_W(s)*species_set.species[s].R;
	T	= p / rhoRmix;



	double E, Er;


	// 1. Enthalpy of formation
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_W(s)*species_set.species[s].h0;


	// 2. Kinetic energy
	double rho_mix	= 0.0;
//#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_W(s);

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= rho_mix*data_W(k)*data_W(k);
	E_k *= 0.5;


	// 3. Translation
	double E_tra = 0.0;
//#pragma omp parallel for reduction(+:E_tra)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_tra	+= data_W(s) * (species_set.species[s].Cv_tra*T);
	}


	// 4. Rotation
	Er = 0.0;
//#pragma omp parallel for reduction(+:Er)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Er	+= data_W(s) * (species_set.species[s].Cv_rot*Tr);
	}


	E = E_tra + E_h0 + E_k + Er;


	data_Q(indexT)	= E;
	data_Q(indexTr)	= Er;

}

void VariableChangeType5::From_Q(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
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

	double T	= (data_Q(indexT) - E_k - E_h0 - data_Q(indexTr)) / rhoCvtra;
	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type5: fromQ => ALL)");
	}
	//Common::ErrorCheckNonNegative<double>(T, "Temperature(Type1)");

	double Tr	= data_Q(indexTr) / rhoCvrot;

	data_V(indexT)	= T;
	data_V(indexTr)	= Tr;

	data_W(indexT)	= rhoR*T;
	data_W(indexTr)	= Tr;



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
	data_MIX(3) = 1.0 + rhoR/(rhoCvtra+rhoCvrot);
	data_MIX(4) = rhoCvtra;
	data_MIX(5) = rhoCvrot;
}


void VariableChangeType5::From_V(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
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

	double T	= data_V(indexT);
	double Tr	= data_V(indexTr);

	double Er	= rhoCvrot * Tr;
	data_Q(indexT)	= rhoCvtra*T + Er + E_h0 + E_k;
	data_Q(indexTr)	= Tr;

	data_W(indexT)	= rhoR*T;
	data_W(indexTr)	= Tr;


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
	data_MIX(3) = 1.0 + rhoR/(rhoCvtra+rhoCvrot);
	data_MIX(4) = rhoCvtra;
	data_MIX(5) = rhoCvrot;
}



void VariableChangeType5::From_W(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
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
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_W(k)*data_W(k);
	E_k = 0.5*rho*E_k;


	// 3. Calculate T
	int indexT	= species_set.NS + ND;
	int indexTr	= species_set.NS + ND + 1;

	double T	= data_W(indexT) / rhoR;
	double Tr	= data_W(indexTr);

	double Er	= rhoCvrot * Tr;
	data_Q(indexT)	= rhoCvtra*T + Er + E_h0 + E_k;
	data_Q(indexTr)	= Er;

	data_V(indexT)	= T;
	data_V(indexTr)	= Tr;


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
	data_MIX(3) = 1.0 + rhoR/(rhoCvtra+rhoCvrot);
	data_MIX(4) = rhoCvtra;
	data_MIX(5) = rhoCvrot;
}



}
}
