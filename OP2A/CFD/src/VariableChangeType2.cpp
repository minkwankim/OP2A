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
//#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_V(s);


	int indexT	= species_set.NS + ND;
	int indexTe	= indexT + 1;



	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_h0	+= data_V(s)*species_set.species[s].h0;
		}
	}


	// 2. E_k (Kinetic energy)
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
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_int	+= data_V(s) * (species_set.species[s].Cv_tr*data_V(indexT));
		}
	}

	// 4. E_e (Electron energy)
	double E_e = 0.0;
//#pragma omp parallel for reduction(+:E_e)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			E_e	= data_V(s) * (species_set.species[s].Cv_tra*data_V(indexTe) + species_set.species[s].h0);
			break;
		}
	}


	data_Q(indexT)	= E_int + E_k + E_h0 + E_e;
	data_Q(indexTe)	= E_e;
}



void VariableChangeType2::Q_to_V(Data::DataStorage& data_Q, CHEM::SpeciesSet& species_set, int ND, Data::DataStorage& data_V, unsigned int CFD_NT)
{
	double rho_mix	= 0.0;
//#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);


	int indexT	= species_set.NS + ND;
	int indexTe	= indexT + 1;


	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_h0	+= data_Q(s)*species_set.species[s].h0;
		}
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
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCvmix	+= data_Q(s) * species_set.species[s].Cv_tr;
		}
	}

	// 4. rhoCve
	double rhoCve = 0.0;
//#pragma omp parallel for reduction(+:rhoCve)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
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
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			E_h0e	= data_Q(s) * species_set.species[s].h0;
			break;
		}
	}




	double T, Te;
	T	= (data_Q(indexT) - data_Q(indexTe) - E_h0 - E_k) / rhoCvmix;
	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type2: Q => V)");
	}

	Te	= (data_Q(indexTe) - E_h0e) / rhoCve;
	if (Te < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Electron Temperature (Type2: Q => V)");
	}

	//Common::ErrorCheckNonNegative<double>(T,  "Temperature(Type2)");
	//Common::ErrorCheckNonNegative<double>(Te, "Electron Temperature(Type2)");

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
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoRmix	+= data_V(s) * species_set.species[s].R;
		}
	}


	double rhoeRe	= 0.0;
//#pragma omp parallel for reduction(+:rhoeRe)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoeRe	= data_V(s) * species_set.species[s].R;
			break;
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
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoRmix	+= data_W(s) * species_set.species[s].R;
		}
	}


	double rhoeRe	= 0.0;
//#pragma omp parallel for reduction(+:rhoeRe)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoeRe	= data_W(s) * species_set.species[s].R;
			break;
		}
	}


	double T	= (data_W(indexT) - data_W(indexTe)) / rhoRmix;
	double Te	= data_W(indexTe) / rhoeRe;

	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type2: W => V)");
	}
	//Common::ErrorCheckNonNegative<double>(T,  "Temperature(Type2: W=>V)");

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
//#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_Q(s);


	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_h0	+= data_Q(s)*species_set.species[s].h0;
		}
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
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCvmix	+= data_Q(s) * species_set.species[s].Cv_tr;
		}
	}

	// 4. rhoCve
	double rhoCve = 0.0;
//#pragma omp parallel for reduction(+:rhoCve)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
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
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			E_h0e	= data_Q(s) * species_set.species[s].h0;
			break;
		}
	}



	double T, Te;
	T	= (data_Q(indexT) - data_Q(indexTe) - E_h0 - E_k) / rhoCvmix;
	Te	= (data_Q(indexTe) - E_h0e) / rhoCve;

	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type2: Q => W)");
	}

	if (Te < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Electron Temperature (Type2: Q => W)");
	}

	//Common::ErrorCheckNonNegative<double>(T,  "Temperature(Type2: Q =>W)");
	//Common::ErrorCheckNonNegative<double>(Te, "Electron Temperature(Type2): Q => W");


	double rhoRmix	= 0.0;
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoRmix	+= data_Q(s) * species_set.species[s].R;
		}
	}


	double rhoeRe	= 0.0;
//#pragma omp parallel for reduction(+:rhoeRe)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoeRe	= data_Q(s) * species_set.species[s].R;
			break;
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
//#pragma omp parallel for reduction(+:rhoRmix)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoRmix	+= data_W(s) * species_set.species[s].R;
		}
	}

	double rhoeRe	= 0.0;
//#pragma omp parallel for reduction(+:rhoeRe)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			rhoeRe	= data_W(s) * species_set.species[s].R;
			break;
		}
	}


	double T	= (data_W(indexT) - data_W(indexTe)) / rhoRmix;
	double Te	= data_W(indexTe) / rhoeRe;

	//Common::ErrorCheckNonNegative<double>(T,  "Temperature(Type2: W=>V)");
	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type2: W => Q)");
	}


	double rho_mix	= 0.0;
//#pragma omp parallel for reduction(+:rho_mix)
	for (int s = 0; s <= species_set.NS-1; s++)	rho_mix	+= data_W(s);



	// 1. E_h0 (Energy by enthalpy for formation) [Except electrons]
	double E_h0 = 0.0;
//#pragma omp parallel for reduction(+:E_h0)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_h0	+= data_W(s)*species_set.species[s].h0;
		}
	}


	// 2. E_k (Kinetic energy)
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
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			E_int	+= data_W(s) * (species_set.species[s].Cv_tr*T);
		}
	}

	// 4. E_e (Electron energy)
	double E_e = 0.0;
//#pragma omp parallel for reduction(+:E_int)
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type == CHEM::SpeciesType::Electron)
		{
			E_e	= data_W(s) * (species_set.species[s].Cv_tra*Te + species_set.species[s].h0);
			break;
		}
	}


	data_Q(indexT)	= E_int + E_h0 + E_k + E_e;
	data_Q(indexTe)	= E_e;
}


void VariableChangeType2::From_Q(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
{
	// 1. Calculate rho
	double rho = 0.0;
//#pragma omp parallel for reduction(+:rho) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rho	+= data_Q(s);


	// 2. Calculate rhoCv and rhoR
	double rhoCv_wo_e	= 0.0;
	double rhoCve 		= 0.0;
	double rhoR_wo_e	= 0.0;
	double rhoRe		= 0.0;
	double E_h0			= 0.0;
	double E_h0e		= 0.0;

//#pragma omp parallel for reduction(+:rhoCv_wo_e, rhoCve, rhoR_wo_e, rhoRe) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCv_wo_e	+= data_Q(s)*species_set.species[s].Cv_tr;
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
	int indexTe	= species_set.NS + ND + 1;

	double T	= (data_Q(indexT) - data_Q(indexTe) - E_k - E_h0) / rhoCv_wo_e;
	if (T < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Temperature (Type2: fromQ => ALL)");
	}
	//Common::ErrorCheckNonNegative<double>(T, "Temperature(Type1)");

	double Te	= (data_Q(indexTe) - E_h0e) / rhoCve;

	data_V(indexT)	= T;
	data_V(indexTe)	= Te;

	double pe		= rhoRe*Te;
	data_W(indexT)	= rhoR_wo_e*T + pe;
	data_W(indexTe)	= pe;


	// Assign Mixture properties
	// Calculate N
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Xs(s)	= data_Q(s) / species_set.species[s].m;
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
	data_MIX(1) = rhoR_wo_e;
	data_MIX(2) = rho / N * AMU_SI;
	data_MIX(3) = 1.0 + rhoR_wo_e/rhoCv_wo_e;
	data_MIX(4) = rhoCv_wo_e;
}



void VariableChangeType2::From_V(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
{
	// 1. Calculate rho
	double rho = 0.0;
//#pragma omp parallel for reduction(+:rho) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rho	+= data_V(s);


	// 2. Calculate rhoCv and rhoR
	double rhoCv_wo_e	= 0.0;
	double rhoCve 		= 0.0;
	double rhoR_wo_e	= 0.0;
	double rhoRe		= 0.0;
	double E_h0			= 0.0;
	double E_h0e		= 0.0;

//#pragma omp parallel for reduction(+:rhoCv_wo_e, rhoCve, rhoR_wo_e, rhoRe) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCv_wo_e	+= data_V(s)*species_set.species[s].Cv_tr;
			rhoR_wo_e	+= data_V(s)*species_set.species[s].R;

			E_h0		+= data_V(s)*species_set.species[s].h0;
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
	int indexTe	= species_set.NS + ND + 1;
	double	T	= data_V(indexT);
	double	Te	= data_V(indexTe);

	double Ee	= rhoCve*Te + E_h0e;
	double E	= rhoCv_wo_e*T + E_k + E_h0 + Ee;


	data_Q(indexT)	= E;
	data_Q(indexTe)	= Ee;

	double	pe	= rhoRe*Te;
	data_W(indexT)	= rhoR_wo_e*T + pe;
	data_W(indexTe)	= pe;


	// Assign Mixture properties
	// Calculate N

#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Xs(s) = data_Q(s) / species_set.species[s].m;
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
	data_MIX(1) = rhoR_wo_e;
	data_MIX(2) = rho / N * AMU_SI;
	data_MIX(3) = 1.0 + rhoR_wo_e/rhoCv_wo_e;
	data_MIX(4) = rhoCv_wo_e;
}


void VariableChangeType2::From_W(Data::DataStorage& data_Q, Data::DataStorage& data_V, Data::DataStorage& data_W, Data::DataStorage& data_MIX, Data::DataStorage& data_Xs, Data::DataStorage& data_Ys, CHEM::SpeciesSet& species_set, int ND, unsigned int CFD_NT)
{
	// 1. Calculate rho
	double rho = 0.0;
//#pragma omp parallel for reduction(+:rho) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)	rho	+= data_W(s);


	// 2. Calculate rhoCv and rhoR
	double rhoCv_wo_e	= 0.0;
	double rhoCve 		= 0.0;
	double rhoR_wo_e	= 0.0;
	double rhoRe		= 0.0;
	double E_h0			= 0.0;
	double E_h0e		= 0.0;

//#pragma omp parallel for reduction(+:rhoCv_wo_e, rhoCve, rhoR_wo_e, rhoRe) num_threads(CFD_NT)
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			rhoCv_wo_e	+= data_W(s)*species_set.species[s].Cv_tr;
			rhoR_wo_e	+= data_W(s)*species_set.species[s].R;

			E_h0		+= data_W(s)*species_set.species[s].h0;
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
	int indexTe	= species_set.NS + ND + 1;
	double	T	= (data_W(indexT) - data_W(indexTe)) / rhoR_wo_e;
	double	Te	= data_W(indexTe) / rhoRe;


	double Ee	= rhoCve*Te + E_h0e;
	double E	= rhoCv_wo_e*T + E_k + E_h0 + Ee;


	data_Q(indexT)	= E;
	data_Q(indexTe)	= Ee;

	data_V(indexT)	= T;
	data_V(indexTe)	= Te;


	// Assign Mixture properties
	// Calculate N
#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		data_Xs(s) = data_Q(s) / species_set.species[s].m;
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
	data_MIX(1) = rhoR_wo_e;
	data_MIX(2) = rho / N * AMU_SI;
	data_MIX(3) = 1.0 + rhoR_wo_e/rhoCv_wo_e;
	data_MIX(4) = rhoCv_wo_e;
}


}
}
