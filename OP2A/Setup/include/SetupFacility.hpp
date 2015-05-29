/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 26, 2015
 *      			Author: Minkwan Kim
 *
 * SetupFacility.hpp
 * 			-  
 *  
 */
#ifndef SETUPFACILITY_HPP_
#define SETUPFACILITY_HPP_



#include <set>

#include "Common/include/NonCopyable.hpp"
#include "Setup/include/SetupArgs.hpp"

namespace OP2A{
namespace Setup{

class SetupObject;


class Setup_API SetupFacility : public Common::NonCopyable<SetupFacility>
{

public: // methods
	/// @return the instance of this singleton
	static SetupFacility& getInstance();

	/// Registers this object in the facility
	/// If object already exists nothing happens
	void registSetupObj( SetupObject* cobj );

	/// Remove sthe object form the facility
	void unregistSetupObj( SetupObject* cobj );

	/// Configure the dynamic options with the arguments passed
	void setupDynamicOptions ( SetupArgs& args );

private:
	/// Default constructor
	SetupFacility();

	/// Default destructor
	~SetupFacility();

private:
	/// storage of the config objects that are registered
	std::set<SetupObject*>	m_config_objs;

};


}
}



#endif /* SETUPFACILITY_HPP_ */
