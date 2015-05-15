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


#ifndef CONFIGURE_CONGIGARGS_HPP
#define CONFIGURE_CONGIGARGS_HPP

#include "Common/include/OP2A.hpp"
#include "Configure/include/Config.hpp"


namespace OP2A{
namespace Config{

typedef std::string ConfigKey;
typedef std::string ConfigValue;
typedef std::map<ConfigKey, ConfigValue> ConfigMap;

/// 	This class represents an Arguments map.
/// 	It is used for storing arguments in Option parsing.
class CONFIG_API	ConfigArgs : public ConfigMap
{
public:
	/// Default constructor without arguments.
	ConfigArgs();

	/// Default destructor.
	~ConfigArgs();


	/// Get these arguments in the form of string
	/// @returns string with the arguments
	std::string str() const;

	/// Assignement operator for a map of strings
	ConfigArgs& operator=(const std::map<std::string, std::string>& value);

	/// Consume the configurations associated with the keys
	void consume (const std::vector<ConfigKey>& keys);

	/// Remove the configuration arguments that start with the provided key
	void remove_filter (const ConfigKey& key);

	/// Only keep the configuration arguments that start with the provided key
	void pass_filter (const ConfigKey& key);
}; // end of class ConfigArgs



}
}
#endif
