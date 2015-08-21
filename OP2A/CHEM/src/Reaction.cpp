/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 16, 2015
 *      			Author: Minkwan Kim
 *
 * Reaction.cpp
 * 			-  
 *  
 */



#include <limits>

#include "Common/include/Exception_DimensionMatch.hpp"
#include "Common/include/Exception_NaNValue.hpp"
#include "Common/include/Exception_NegativeValue.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"
#include "Common/include/StringOps.hpp"

#include "Math/include/OP2A_Math.hpp"

#include "CHEM/include/ChemConstants.hpp"
#include "CHEM/include/Reaction.hpp"



using namespace std;

namespace OP2A{
namespace CHEM{


double Calc_Tc(double T, double Tr, double Tv, double Te, ReactionType type, int direction)
{
	double Tc;
	double Ttr = pow(T, 1/2.0) * pow(Tr, 1/2.0);
	double Td  = pow(T, 1/4.0) * pow(Tr, 1/4.0) * pow(Tv, 1/4.0) * pow(Te, 1/4.0);;

	if(direction == 0)
	{
		switch (type)
		{
		case ReactionType::DISSOCIATION:
			Tc = Td;
			break;

		case ReactionType::EXCHANGE:
			Tc = Ttr;
			break;

		case ReactionType::DISSOCIATIVE_RECOMBINATION:
			Tc = Ttr;
			break;

		case ReactionType::CHARGE_EXCHANGE:
			Tc = Ttr;
			break;

		case ReactionType::ELECTRON_IMPACT_DISSOCIATION:
			Tc = Te;
			break;

		case ReactionType::ELECTRON_IMPACT_IONIZATION:
			Tc = Te;
			break;
		}
	}
	else
	{
		switch (type)
		{
		case ReactionType::DISSOCIATION:
			Tc = Td;
			break;

		case ReactionType::EXCHANGE:
			Tc = Ttr;
			break;

		case ReactionType::DISSOCIATIVE_RECOMBINATION:
			Tc = Te;
			break;

		case ReactionType::CHARGE_EXCHANGE:
			Tc = Ttr;
			break;

		case ReactionType::ELECTRON_IMPACT_DISSOCIATION:
			Tc = Te;
			break;

		case ReactionType::ELECTRON_IMPACT_IONIZATION:
			Tc = Te;
			break;
		}
	}

	return (Tc);
}



Reaction::Reaction():ID(-1), type(-1), m_allocated(false), m_hasdata(false), method_kb(0), Tref(0.0)
{

}

Reaction::Reaction(int NS):ID(-1), type(-1), m_allocated(false), m_hasdata(false), method_kb(0), Tref(0.0)
{
	alpha.resize(NS, 0.0);
	beta.resize(NS, 0.0);
	beta_m_alpha.resize(NS, 0.0);
}


Reaction::Reaction(std::string data_react, std::string data_prod):ID(-1), type(-1), m_allocated(true), m_data_react(data_react), m_data_prod(data_prod), m_hasdata(true), method_kb(0), Tref(0.0)
{

}


Reaction::Reaction(std::string data_react, std::string data_prod, int i_type):ID(-1), type(i_type), m_allocated(true), m_data_react(data_react), m_data_prod(data_prod), m_hasdata(true),  method_kb(0), Tref(0.0)
{

}



Reaction::Reaction(std::string data_raw):ID(-1), type(-1), m_allocated(false)
{
	data_assigning(data_raw);
	data_processing();
}


Reaction::Reaction(std::string data_raw, int i_type):ID(-1), type(i_type), m_allocated(false)
{
	data_assigning(data_raw);
	data_processing();
}






Reaction::~Reaction()
{

}



/*
 * Internal function
 */
void Reaction::data_processing()
{
	if (m_hasdata == false)
	{
		throw Common::ExceptionDimensionMatch (FromHere(), "Problem in Reaction DATA: Need raw data to process it");
	}

	Common::StringOps::split(m_data_react, '+', m_React);
	Common::StringOps::split(m_data_prod, '+', m_Prod);

	m_React_c.resize(m_React.size());
	m_Prod_c.resize(m_Prod.size());

	for (int i = 0; i <= m_React.size()-1; i++)
	{
		Common::StringOps::trim(m_React[i]);
		m_React_c[i] = Common::StringOps::extractingNumber<double>(m_React[i]);
	}

	for (int i = 0; i <= m_Prod.size()-1; i++)
	{
		Common::StringOps::trim(m_Prod[i]);
		m_Prod_c[i] = Common::StringOps::extractingNumber<double>(m_Prod[i]);
	}
}

void Reaction::data_processing(std::string data_raw)
{
	data_assigning (data_raw);

	Common::StringOps::split(m_data_react, '+', m_React);
	Common::StringOps::split(m_data_prod, '+', m_Prod);

	m_React_c.resize(m_React.size());
	m_Prod_c.resize(m_Prod.size());

	for (int i = 0; i <= m_React.size()-1; i++)
	{
		Common::StringOps::trim(m_React[i]);
		m_React_c[i] = Common::StringOps::extractingNumber<double>(m_React[i]);
	}

	for (int i = 0; i <= m_Prod.size()-1; i++)
	{
		Common::StringOps::trim(m_Prod[i]);
		m_Prod_c[i] = Common::StringOps::extractingNumber<double>(m_Prod[i]);
	}
}



void Reaction::data_assigning(std::string data_raw)
{
	std::vector<std::string> str_temp;
	Common::StringOps::split(data_raw, '=', str_temp);

	if (str_temp.size() != 2)
	{
		throw Common::ExceptionDimensionMatch (FromHere(), "Problem in Reaction DATA: Reactant/Product separation");
	}

	Common::StringOps::trim(str_temp[0]);
	Common::StringOps::trim(str_temp[1]);

	m_data_react = str_temp[0];
	m_data_prod  = str_temp[1];
	m_hasdata = true;
}


void Reaction::data_completing(const std::vector<Species>& species)
{
	int NS	= species.size();
	alpha.resize(NS, 0.0);
	beta.resize(NS, 0.0);
	beta_m_alpha.resize(NS, 0.0);

	for (int s = 0; s <= NS-1; s++)
	{
		alpha[s]		= 0.0;
		beta[s]			= 0.0;
		beta_m_alpha[s] = 0.0;
	}

	for (int s1 = 0; s1 <= m_React.size()-1; s1++)
	{
		for (int s2 = 0; s2 <= NS-1; s2++)
		{
			if (m_React[s1] == species[s2].name)
			{
				alpha[s2]	+= m_React_c[s1];
			}
		}
	}

	for (int s1 = 0; s1 <= m_Prod.size()-1; s1++)
	{
		for (int s2 = 0; s2 <= NS-1; s2++)
		{
			if (m_Prod[s1] == species[s2].name)
			{
				beta[s2]	+= m_Prod_c[s1];
			}
		}
	}

	for (int s = 0; s <= NS-1; s++)
	{
		beta_m_alpha[s] = beta[s] - alpha[s];
	}

	m_allocated = true;
}


double Reaction::kf(double T)
{
	return (ForwardReaction.ReactionRate(T));
}

double Reaction::kf(double T, double n_mix)
{
	return (ForwardReaction.ReactionRate(T));
}


double Reaction::kb(double T)
{
	return (BackwardReaction.ReactionRate(T));
}

double Reaction::kb(double T, double n_mix)
{
	if (method_kb == 0)
	{
		return (kb(T));
	}
	else
	{
		double keq = Keq.Keq(T, n_mix);
		double kfb = kf(T);
		return (kfb/keq);
	}
}




}
}
