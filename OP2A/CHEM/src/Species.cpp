/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 11, 2015
 *      			Author: Minkwan Kim
 *
 * Species.cpp
 * 			-  
 *  
 */


#include <limits>

#include "Common/include/Exception_DimensionMatch.hpp"
#include "Common/include/Exception_NaNValue.hpp"
#include "Common/include/Exception_NegativeValue.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"

#include "CHEM/include/ChemConstants.hpp"
#include "CHEM/include/Species.hpp"
#include "Math/include/OP2A_Math.hpp"



using namespace std;

namespace OP2A{
namespace CHEM{



Species::Species(): data_assigned(false)
{

}

Species::Species(const std::string& species_name)
{
	AssignData(species_name);
}

Species::~Species()
{

}




void Species::AssignData(const std::string& species_name)
{
	AssignDataBasic(species_name);
	AssignThermoData(species_name, R, type);
	AssignTransport(species_name);
	data_assigned = true;
}

bool Species::is_data_assigned()
{
	return (data_assigned);
}








/*
 * Internal functions
 */
double Species::Cv_TRA(const double& T)
{
	return Cv_tra;
}

double Species::Cv_ROT(const double& T)
{
	return Cv_rot;
}

double Species::Cv_TR(const double& T)
{
	return Cv_tr;
}



double Species::Cv_VIB(const double& T)
{
	int n;
	double Cv;
	double aux2, aux3;

	if (T >= 10.0)
	{
		Cv	= 0.0;
		if (type == SpeciesType::Molecule)
		{
			for (n = 0;  n <= n_vib_lvl-1; n++)
			{
				aux2	= theta_vib[n]/T;
				aux3	= exp(aux2);
				Cv		+= R * aux2*aux2*aux3 / pow((aux3-1.0), 2.0);
			}
		}

		if (Cv < 0.0)
		{
			throw Common::ExceptionNegativeValue (FromHere(), "Negative Value for Cv_VIB");
		}

		if (Cv == numeric_limits<double>::infinity())
		{
			throw Common::ExceptionInfiniteValue (FromHere(), "Inifnite Value for Cv_VIB");
		}

		if (Cv != Cv)
		{
			throw Common::ExceptionNaNValue (FromHere(), "Inifnite Value for Cv_VIB");
		}
	}
	else
	{
		Cv = 0.0;
	}

	return (Cv);
}


double Species::Cv_ELE(const double& T)
{
	int n;
	double Cv;
	double aux2, aux3, aux4;
	double component1, component2, component3, component4;

	component1	= 0.0;
	component2	= 0.0;
	component3	= 0.0;
	component4	= 0.0;

	Cv	= 0.0;
	if (type != SpeciesType::Electron)
	{
		aux2		= theta_el[0]/ T;
		aux3		= exp(-aux2);

		component1	+= g_el[0] * aux3;
		component2 	+= g_el[0] * aux2/T * aux3;

		for (n = 1;  n <= n_elec_lvl; n++)
		{
			aux2		= theta_el[n]/ T;
			aux3		= exp(-aux2);

			component1	+= g_el[n] * aux3;
			component2 	+= g_el[n] * aux2/T * aux3;
			component3	+= g_el[n] * aux2*aux2 * aux3;
			component4  += g_el[n] * theta_el[n] * aux3;
		}

		Cv	= R * ((component3/component1) - (component4*component2)/(component1*component1));
	}

	if (Cv < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Value for Cv_ELE");
	}

	if (Cv == numeric_limits<double>::infinity())
	{
		throw Common::ExceptionInfiniteValue (FromHere(), "Inifnite Value for Cv_ELE");
	}

	if (Cv != Cv)
	{
		throw Common::ExceptionNaNValue (FromHere(), "Inifnite Value for Cv_ELE");
	}

	return (Cv);
}

double Species::Cv_VE(const double& T)
{
	double Cv;

	Cv = Cv_VIB(T)	+ Cv_ELE(T);
	return (Cv);
}


double Species::Cv_VEE(const double& T)
{
	double Cv;

	if (type != SpeciesType::Electron)	Cv = Cv_VE(T);
	else								Cv = Cv_tra;

	return (Cv);
}



double Species::e_VIB(const double& T)
{
	int n;
	double evib = 0.0;
	if (type == SpeciesType::Molecule)
	{
		if (T >= 10.0)
		{
			for (n = 0;  n <= n_vib_lvl-1; n++)
			{
				evib	+= R * theta_vib[n] / (exp(theta_vib[n]/T) - 1.0);
			}
		}
	}

	if (evib < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Value for e_vib");
	}

	if (evib == numeric_limits<double>::infinity())
	{
		throw Common::ExceptionInfiniteValue (FromHere(), "Inifnite Value for e_vib");
	}

	if (evib != evib)
	{
		throw Common::ExceptionNaNValue (FromHere(), "NaN Value for e_vib");
	}


	return(evib);
}

double Species::e_ELE(const double& T)
{
	int n;
	double e_el;
	double aux2, aux3, aux4;

	e_el	= 0.0;
	if (type != SpeciesType::Electron)
	{
		aux2	= 0.0;
		aux3	= 0.0;

		if (n_elec_lvl > 0)
		{
			for (n = 0;  n <= n_elec_lvl-1; n++)
			{
				aux4	= g_el[n]	* exp(-theta_el[n]/ T);
				aux2	+= aux4 * theta_el[n];
				aux3	+= aux4;
			}

			e_el	= R * aux2/aux3;
		}
	}


	if (e_el < 0.0)
	{
		throw Common::ExceptionNegativeValue (FromHere(), "Negative Value for e_ELE");
	}

	if (e_el == numeric_limits<double>::infinity())
	{
		throw Common::ExceptionInfiniteValue (FromHere(), "Inifnite Value for e_ELE");
	}

	if (e_el != e_el)
	{
		throw Common::ExceptionNaNValue (FromHere(), "NaN Value for e_ELE");
	}

	return (e_el);
}

double Species::e_VE(const double& T)
{
	double eve;

	eve = e_VIB(T)	+ e_ELE(T);
	return (eve);
}


double Species::e_VEE(const double& T)
{
	double evee;

	if (type != SpeciesType::Electron)	evee = e_VE(T);
	else								evee = Cv_VEE(T) * T;

	return (evee);
}



double Species::e_internal(const double& Ttra, const double& Trot, const double& Tvib, const double& Tele)
{
	double es;

	es	= Cv_tra*Ttra;
	es	+= Cv_rot*Trot;
	es	+= e_VIB(Tvib);
	es	+= e_ELE(Tele);

	return (es);
}


double Species::e_internal(const double& T)
{
	double es;

	es	= Cv_tra*T;
	es	+= Cv_rot*T;
	es	+= e_VIB(T);
	es	+= e_ELE(T);

	return (es);
}

double Species::e_internal(const double& Ttra, const double& Tvib)
{
	double es;

	es	= Cv_tra*Ttra;
	es	+= Cv_rot*Ttra;
	es	+= e_VIB(Tvib);
	es	+= e_ELE(Tvib);

	return (es);
}

double Species::e_internal(const double& Ttra, const double& Trot, const double& Tvib)
{
	double es;

	es	= Cv_tra*Ttra;
	es	+= Cv_rot*Trot;
	es	+= e_VIB(Tvib);
	es	+= e_ELE(Tvib);

	return (es);
}


double Species::h(const double& Ttra, const double& Trot, const double& Tvib, const double& Tele)
{
	double h;
	double es	= e_internal(Ttra, Trot, Tvib, Tele);

	h	= es + R*Ttra;

	return (h);
}

double Species::h(const double& Ttra)
{
	double h;
	double es	= e_internal(Ttra);

	h	= es + R*Ttra;

	return (h);
}


double Species::h(const double& Ttra, const double& Tvib)
{
	double h;
	double es	= e_internal(Ttra, Tvib);

	h	= es + R*Ttra;

	return (h);
}


double Species::h(const double& Ttra, const double& Trot, const double& Tvib)
{
	double h;
	double es	= e_internal(Ttra, Trot, Tvib);

	h	= es + R*Ttra;

	return (h);
}




double Species::viscosityCoeff(const double& T, unsigned int & model)
{
	double mu;
	double temp1, temp2, temp3;

	if (type == SpeciesType::Electron)
	{
		mu	= 0.0;
	}
	else
	{
		if (model == SpeciesViscosityModel::VHS_MODEL)
		{
			temp1	= sqrt(MATH_PI * C_BOLTZMANN_SI * m * T_ref);
			temp2	= 2.0 * MATH_PI * pow(d_ref, 2.0) * (5.0-2.0*omega) * (7.0-2.0*omega);
			mu = 15.0 * temp1 / temp2 * pow(T/T_ref, omega);
		}
		else if (model	== SpeciesViscosityModel::BLOTTNER_MODEL)
		{
			double lnT = log(T);
			temp1	= Blottner[0]*lnT + Blottner[1];
			mu 		= 0.1 * exp(temp1*lnT + Blottner[2]);
		}
		else if (model	== 	SpeciesViscosityModel::SUTHERLAND_MODEL)
		{
			temp1 	= Sutherland[1] + Sutherland[2];
			temp2 	= T + Sutherland[2];
			temp3 	= T / Sutherland[1];

			mu		= Sutherland[0] * (temp1/temp2) * pow(temp3, 1.5);
		}
		else if (model	== SpeciesViscosityModel::KINETIC_MODEL)
		{
			temp1	= 1.147 * pow(T/T_eps, -0.145) + pow(T/T_eps + 0.5, -2.0);
			mu		= (2.68E-6 * sqrt(M*T)) / (pow(sigma, 2.0)*temp1);
		}
	}


	if (mu	!= mu) throw Common::ExceptionNaNValue (FromHere(), "Nan value for Viscosity coefficient ");
	if (mu == numeric_limits<double>::infinity()) throw Common::ExceptionInfiniteValue (FromHere(), "Inifnite Value for viscosity coefficient");


	return (mu);

}



double  Species::thermalConductivityCoeff(const double& Tvib, const double& Tele, const double& mu,int model, int mode)
{
	double kappa;

	switch(model)
	{
	case SpeciesThermalConductivityModel::EUCKEN_MODEL:
		if (mode == EnergyMode::E_TRA)
		{
			kappa	= 2.5 * mu * Cv_tra;
		}
		else if (mode == EnergyMode::E_ROT)
		{
			kappa	= mu * Cv_rot;
		}
		else if (mode == EnergyMode::E_VIB)
		{
			kappa	= mu * Cv_VIB(Tvib);
		}
		else if (mode == EnergyMode::E_ELE)
		{
			kappa	= mu * Cv_VIB(Tele);
		}
		else if (mode == EnergyMode::E_TR)
		{	double kappa_tra	= thermalConductivityCoeff(Tvib, Tele, mu, model, 0);
			double kappa_rot	= thermalConductivityCoeff(Tvib, Tele, mu, model, 1);
			kappa	= kappa_tra + kappa_rot;
		}
		else if (mode == EnergyMode::E_VE)
		{
			double kappa_vib	= thermalConductivityCoeff(Tvib, Tele, mu, model, 2);
			double kappa_ele	= thermalConductivityCoeff(Tvib, Tvib, mu, model, 3);
			kappa	= kappa_vib + kappa_ele;
		}
		else if (mode == EnergyMode::E_VEE)
		{
			if (type != SpeciesType::Electron)	kappa	= thermalConductivityCoeff(Tvib, Tele, mu, model, 5);
			else								kappa	= thermalConductivityCoeff(Tvib, Tele, mu, model, 0) / 2.5;
		}
		break;
	}


	if (kappa	!= kappa) throw Common::ExceptionNaNValue (FromHere(), "Nan value for Thermal Conductivity ");
	if (kappa == numeric_limits<double>::infinity()) throw Common::ExceptionInfiniteValue (FromHere(), "Inifnite Value for Thermal Conductivity ");


	return (mu);

}





}
}
