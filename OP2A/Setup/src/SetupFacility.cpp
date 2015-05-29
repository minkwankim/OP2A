/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 26, 2015
 *      			Author: Minkwan Kim
 *
 * SetupFacility.cpp
 * 			-  
 *  
 */


#include "Setup/include/SetupFacility.hpp"
#include "Setup/include/SetupObject.hpp"


using namespace std;

namespace OP2A{
namespace Setup{


/*
 * Constructor and destructor
 */
SetupFacility::SetupFacility()
{

}

SetupFacility::~SetupFacility()
{

}



SetupFacility& SetupFacility::getInstance()
{
	static SetupFacility theSetupFacility;
	return theSetupFacility;
}

void SetupFacility::registSetupObj(SetupObject* cobj)
{
	op_assert( cobj != OPNULL );
	m_config_objs.insert(cobj);
}

void SetupFacility::unregistSetupObj(SetupObject* cobj)
{
	op_assert( cobj != OPNULL );
	m_config_objs.erase(cobj);
}

void SetupFacility::setupDynamicOptions(SetupArgs& args)
{
	for ( std::set<SetupObject*>::iterator itr = m_config_objs.begin(), stop = m_config_objs.end(); itr != stop; ++itr )
	{
		SetupObject * cobj = *itr;
		op_assert ( *itr != OPNULL );
		cobj->setupDynamic(args);
	}
}

}
}


