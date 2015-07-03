/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 16, 2015
 *      			Author: Minkwan Kim
 *
 * VariableChange_Type1.cpp
 * 			-  
 *  
 */


#include <omp.h>

#include "CFD/include/VariableChangeTypes.hpp"


namespace OP2A{
namespace CFD{



void VariableChangeType1::V_to_Q(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	double rho_mix	= 0.0;
//#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_V(s);


	// Total Energy
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_V(s)*species_set.species[s].h0;

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= rho_mix*data_V(k)*data_V(k);
	E_k *= 0.5;


	int indexT	= species_set.NS + ND;
	double E_int = 0.0;
//#pragma omp parallel for reduction(+:E_int)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_int	+= data_V(s) * (species_set.species[s].Cv_tr*data_V(indexT));
		//E_int	+= data_V(s) * (species_set.species[s].Cv_tr*data_V(indexT) + species_set.species[s].e_VIB(data_V(indexT)) + species_set.species[s].e_ELE(data_V(indexT)));
	}

	data_Q(indexT)	= E_int + E_k + E_h0;
}



void VariableChangeType1::Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	double rho_mix	= 0.0;
//#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);


	// Translational Temperature
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_Q(s)*species_set.species[s].h0;

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_Q(k)*data_Q(k) / rho_mix;
	E_k *= 0.5;


	int indexT	= species_set.NS + ND;

	double Cv_mix = 0.0;
//#pragma omp parallel for reduction(+:Cv_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Cv_mix	+= data_Q(s) * species_set.species[s].Cv_tr;
	}

	double T	= (data_Q(indexT) - E_k - E_h0) / Cv_mix;

	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type1: Q => V)");
	}
	//Common::ErrorCheckNonNegative<double>(T, "Temperature(Type1)");


	data_V(indexT)	= T;
}



//////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * V ><=> W
 */

void VariableChangeType1::V_to_W(Data::DataStorage& data_V, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_V(s)*species_set.species[s].R;


	// Pressure
	int indexT	= species_set.NS + ND;
	data_W(indexT)	= rhoRmix*data_V(indexT);
}


void VariableChangeType1::W_to_V(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_W(s)*species_set.species[s].R;


	// Temperature
	int indexT	= species_set.NS + ND;
	data_V(indexT)	= data_W(indexT) / rhoRmix;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Q <=> W
 */
void VariableChangeType1::Q_to_W(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_W, unsigned int CFD_NT)
{
	double rho_mix	= 0.0;
//#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);


	// Translational Temperature
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_Q(s)*species_set.species[s].h0;

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_Q(k)*data_Q(k) / rho_mix;
	E_k *= 0.5;


	int indexT	= species_set.NS + ND;
	double Cv_mix = 0.0;
//#pragma omp parallel for reduction(+:Cv_mix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Cv_mix	+= data_Q(s) * species_set.species[s].Cv_tr;
	}

	double T	= (data_Q(indexT) - E_k - E_h0) / Cv_mix;


	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_Q(s)*species_set.species[s].R;


	// Pressure
	data_W(indexT)	= rhoRmix*T;
}


void VariableChangeType1::W_to_Q(Data::DataStorage& data_W, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_Q, unsigned int CFD_NT)
{
	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoRmix	+= data_W(s)*species_set.species[s].R;

	// Temperature
	int indexT	= species_set.NS + ND;
	double T = data_W(indexT) / rhoRmix;

	double rho_mix	= 0.0;
//#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_W(s);


	// Total Energy
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_W(s)*species_set.species[s].h0;

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= rho_mix*data_W(k)*data_W(k);
	E_k *= 0.5;

	double E_int = 0.0;
//#pragma omp parallel for reduction(+:E_int)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		E_int	+= data_W(s) * (species_set.species[s].Cv_tr*T);
	}

	data_Q(indexT)	= E_int + E_k + E_h0;
}



void VariableChangeType1::From_Q(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
{
	// 1. Calculate rho
	double rho = 0.0;
//#pragma omp parallel for reduction(+:rho)
	for (int s = 0; s <= species_set.NS-1; s++)	rho	+= data_Q.data[s];


	// 2. Calculate rhoCv and rhoR
	double rhoCv = 0.0;
	double rhoR  = 0.0;

//#pragma omp parallel for reduction(+:rhoCv)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoCv	+= data_Q.data[s]*species_set.species[s].Cv_tr;

//#pragma omp parallel for reduction(+:rhoR)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoR	+= data_Q.data[s]*species_set.species[s].R;



	// 1. Species and velocities
//#pragma omp parallel for num_threads(CFD_NT)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_V.data[s]	= data_Q.data[s];
		data_W.data[s]	= data_Q.data[s];
	}

//#pragma omp parallel for num_threads(CFD_NT)
#pragma ivdep
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		data_V.data[k]	= data_Q.data[k] / rho;
		data_W.data[k]	= data_V.data[k];
	}

	// 2. Calculate E_ho and E_k
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_Q.data[s]*species_set.species[s].h0;

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_V.data[k]*data_V.data[k];
	E_k = 0.5*rho*E_k;


	// 3. Calculate T
	int indexT	= species_set.NS + ND;
	double T	= (data_Q.data[indexT] - E_k - E_h0) / rhoCv;
	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type1: fromQ => ALL)");
	}
	//Common::ErrorCheckNonNegative<double>(T, "Temperature(Type1)");

	data_V.data[indexT]	= T;
	data_W.data[indexT]	= rhoR*T;



	// Assign Mixture properties
	// Calculate N
	double N = 0.0;

#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Xs.data[s]	=	data_Q.data[s] / species_set.species[s].m;
	}

//#pragma omp parallel for reduction(+:N)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		N	+= data_Xs.data[s];
	}

//#pragma omp parallel for num_threads(CFD_NT)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Xs.data[s]	/= N;
		data_Ys.data[s]	= data_Q.data[s] / rho;
	}


	data_MIX.data[0] = rho;
	data_MIX.data[1] = rhoR;
	data_MIX.data[2] = rho / N * AMU_SI;
	data_MIX.data[3] = 1.0 + rhoR/rhoCv;
	data_MIX.data[4] = rhoCv;
}


void VariableChangeType1::From_V(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
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
	int indexT	= species_set.NS + ND;
	double T	= data_V(indexT);

	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_V(s)*species_set.species[s].h0;

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_V(k)*data_V(k);
	E_k = 0.5*rho*E_k;

	// 3. Calculate E and p
	double E	= rhoCv*T + E_h0 + E_k;
	data_Q(indexT)	= E;
	data_W(indexT)	= rhoR*T;



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
	data_MIX(3) = 1.0 + rhoR/rhoCv;
	data_MIX(4) = rhoCv;
}


void VariableChangeType1::From_W(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
{
	// 1. Calculate rho
	double rho = 0.0;
//#pragma omp parallel for reduction(+:rho)
	for (int s = 0; s <= species_set.NS-1; s++)	rho	+= data_W(s);


	// 2. Calculate rhoCv and rhoR
	double rhoCv = 0.0;
	double rhoR  = 0.0;

//#pragma omp parallel for reduction(+:rhoCv)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoCv	+= data_W(s)*species_set.species[s].Cv_tr;

//#pragma omp parallel for reduction(+:rhoR)
	for (int s = 0; s <= species_set.NS-1; s++)	rhoR	+= data_W(s)*species_set.species[s].R;



	// 1. Species and velocities
#pragma ivdep
//#pragma omp parallel for num_threads(CFD_NT)
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
	int indexT	= species_set.NS + ND;
	double T	= data_W(indexT) / rhoR;

	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)	E_h0	+= data_V(s)*species_set.species[s].h0;

	double E_k = 0.0;
//#pragma omp parallel for reduction(+:E_k)
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)	E_k	+= data_V(k)*data_V(k);
	E_k = 0.5*rho*E_k;

	// 3. Calculate E and p
	double E	= rhoCv*T + E_h0 + E_k;
	data_Q(indexT)	= E;
	data_V(indexT)	= T;



	// Assign Mixture properties
	// Calculate N

#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Xs.data[s]	=	data_Q.data[s] / species_set.species[s].m;
	}

	double N = 0.0;
//#pragma omp parallel for reduction(+:N)
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
