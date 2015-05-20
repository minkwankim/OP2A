/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 18, 2015
 *      			Author: Minkwan Kim
 *
 * SetupArgs.hpp
 * 			-  
 *  
 */
#ifndef SETUPARGS_HPP_
#define SETUPARGS_HPP_


#include "Common/include/OP2A.hpp"
#include "Setup/include/SetupAPI.hpp"


namespace OP2A {
namespace Setup {



/*
 * Type definitions for Setup
 */

// Definition of the SetupKey type
typedef std::string SetupKey;

// Definition of the SetupValue type

typedef std::string SetupValue;

/// Definition of the map std::string to std::string type.
typedef std::map<SetupKey, SetupValue> SetupMap;


/*
 * Class for an Arguments map
 * @author	Minkwan Kim
 * @version 1.0 18/5/2015
 *
 */
class Setup_API SetupArgs : public SetupMap
{
public:

	// Constructor and destructor
	SetupArgs()	{	};
	~SetupArgs() {	};


	//	F01- Get arguments in the form of string
	/// @returns string with the arguments
	std::string str() const;


	// F02- Assignement operator for a map of strings
	SetupArgs& operator=(const std::map<std::string, std::string>& value);

	// F03- Consume the configurations associated with the keys
	void consume (const std::vector<SetupKey>& keys);

	// F04- Remove the configuration arguments that start with the provided key
	void remove_filter (const SetupKey& key);

	// F05- Only keep the configuration arguments that start with the provided key
	void pass_filter (const SetupKey& key);
};


}
}



#endif /* SETUPARGS_HPP_ */
