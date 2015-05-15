/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 11, 2015
 *      			Author: Minkwan Kim

 *	[This is adapted from CoolFluid of VKI] 			-
 *
 */

#include <sstream>

#include "Common/include/OP2A.hpp"
#include "Common/include/StringOps.hpp"

#include "Configure/include/Config.hpp"
#include "Configure/include/ConfigArgs.hpp"



using namespace std;
using namespace OP2A::Common;



namespace OP2A {
namespace Config {




// Constructor and Destructor
ConfigArgs::ConfigArgs()	{ }
ConfigArgs::~ConfigArgs()	{ }



//////////////////////////////////////////////////////////////////////////////
std::string ConfigArgs::str () const
{
	std::ostringstream oss;
	for (ConfigArgs::const_iterator it = begin(); it != end(); ++it)	oss << it->first << " " << it->second << "\n";
	return oss.str();
}


ConfigArgs& ConfigArgs::operator=(const std::map<std::string, std::string>& value)
{
	std::map<std::string, std::string>::operator=(value);
	return *this;
}


void ConfigArgs::consume(const std::vector<ConfigKey>& keys)
{
	std::vector<ConfigKey>::const_iterator itr = keys.begin();
	for (; itr != keys.end(); ++itr)	erase( *itr );
}


void ConfigArgs::remove_filter(const ConfigKey& filtkey)
{
	std::vector<ConfigKey> remove_args;
	for (ConfigArgs::iterator it = begin(); it != end(); ++it)
	{
		const ConfigKey& this_key = it->first;

		if (StringOps::startsWith(this_key, filtkey))	remove_args.push_back(this_key);
	}
	consume(remove_args);
}


void ConfigArgs::pass_filter(const ConfigKey& filtkey)
{
	std::vector<ConfigKey> remove_args;
	for (ConfigArgs::iterator it = begin(); it != end(); ++it)
	{
		const ConfigKey& this_key = it->first;

		if ( ! StringOps::startsWith(this_key,filtkey) )	remove_args.push_back(this_key);
	}
	consume(remove_args);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

