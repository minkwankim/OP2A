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

Reaction::Reaction():ID(-1), type(-1), m_allocated(false), m_hasdata(false)
{

}

Reaction::Reaction(int NS):ID(-1), type(-1), m_allocated(false), m_hasdata(false)
{
	alpha.resize(NS, 0.0);
	beta.resize(NS, 0.0);
	alpha_m_beta.resize(NS, 0.0);
}


Reaction::Reaction(std::string data_react, std::string data_prod):ID(-1), type(-1), m_allocated(true), m_data_react(data_react), m_data_prod(data_prod), m_hasdata(true)
{

}


Reaction::Reaction(std::string data_react, std::string data_prod, int i_type):ID(-1), type(i_type), m_allocated(true), m_data_react(data_react), m_data_prod(data_prod), m_hasdata(true)
{

}



Reaction::Reaction(std::string data_raw):ID(-1), type(-1), m_allocated(false)
{
	std::vector<std::string> str_temp;
	Common::StringOps::split(data_raw, '=', str_temp);

	if (str_temp.size() != 2)
	{
		throw Common::ExceptionDimensionMatch (FromHere(), "Problem in Reaction DATA: Reactant/Product separation");
	}

	m_data_react = str_temp[0];
	m_data_prod  = str_temp[1];
	m_hasdata = true;

	data_processing();
}


Reaction::Reaction(std::string data_raw, int i_type):ID(-1), type(i_type), m_allocated(false)
{
	std::vector<std::string> str_temp;
	Common::StringOps::split(data_raw, '=', str_temp);

	if (str_temp.size() != 2)
	{
		throw Common::ExceptionDimensionMatch (FromHere(), "Problem in Reaction DATA: Reactant/Product separation");
	}

	m_data_react = str_temp[0];
	m_data_prod  = str_temp[1];
	m_hasdata = true;

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
		m_React_c[i] = Common::StringOps::extractingNumber<double>(m_React[i]);
	}

	for (int i = 0; i <= m_Prod.size()-1; i++)
	{
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
		m_React_c[i] = Common::StringOps::extractingNumber<double>(m_React[i]);
	}

	for (int i = 0; i <= m_Prod.size()-1; i++)
	{
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

	m_data_react = str_temp[0];
	m_data_prod  = str_temp[1];
	m_hasdata = true;
}


void Reaction::data_completing(const std::vector<Species>& species)
{
	int NS	= species.size();
	alpha.resize(NS, 0.0);
	beta.resize(NS, 0.0);
	alpha_m_beta.resize(NS, 0.0);

	for (int s = 0; s <= NS-1; s++)
	{
		alpha[s]		= 0.0;
		beta[s]			= 0.0;
		alpha_m_beta[s] = 0.0;
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
		alpha_m_beta[s] = alpha[s] - beta[s];
	}

}


}
}
