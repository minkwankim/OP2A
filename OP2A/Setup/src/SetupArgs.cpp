/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 18, 2015
 *      			Author: Minkwan Kim
 *
 * SetupArgs.cpp
 * 			-  
 *  
 */



#include <sstream>

#include "Common/include/OP2A.hpp"
#include "Common/include/StringOps.hpp"

#include "Setup/include/SetupAPI.hpp"
#include "Setup/include/SetupArgs.hpp"


using namespace std;
using namespace OP2A::Common;


namespace OP2A {
namespace Setup {


// F01- Get arguments in the form of string
std::string SetupArgs::str () const
{
	std::ostringstream oss;
	for (SetupArgs::const_iterator it = begin(); it != end(); ++it)	oss << it->first << " " << it->second << "\n";
	return oss.str();
}


// F02- Assignement operator for a map of strings
SetupArgs& SetupArgs::operator=(const std::map<std::string, std::string>& value)
{
	std::map<std::string, std::string>::operator=(value);
	return *this;
}

// F03- Consume the configurations associated with the keys
void SetupArgs::consume(const std::vector<SetupKey>& keys)
{
	std::vector<SetupKey>::const_iterator itr = keys.begin();
	for (; itr != keys.end(); ++itr)	erase( *itr );
}


// F04- Remove the configuration arguments that start with the provided key
void SetupArgs::remove_filter(const SetupKey& filtkey)
{
	std::vector<SetupKey> remove_args;
	for (SetupArgs::iterator it = begin(); it != end(); ++it)
	{
		const SetupKey& this_key = it->first;

		if (StringOps::startsWith(this_key, filtkey))	remove_args.push_back(this_key);
	}
	consume(remove_args);
}


// F05- Only keep the configuration arguments that start with the provided key
void SetupArgs::pass_filter(const SetupKey& filtkey)
{
	std::vector<SetupKey> remove_args;
	for (SetupArgs::iterator it = begin(); it != end(); ++it)
	{
		const SetupKey& this_key = it->first;

		if ( ! StringOps::startsWith(this_key,filtkey) )	remove_args.push_back(this_key);
	}
	consume(remove_args);
}


}
}

