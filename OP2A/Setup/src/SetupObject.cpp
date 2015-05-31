/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 26, 2015
 *      			Author: Minkwan Kim
 *
 * SetupObject.cpp
 * 			-  
 *  
 */

#include <cstdlib>
#include <fstream>

#include "Common/include/FileUtilities.hpp"
#include "Common/include/Exception_FileSystem.hpp"

#include "Setup/include/SetupFileReader.hpp"
#include "Setup/include/SetupFacility.hpp"
#include "Setup/include/SetupObject.hpp"
#include "Setup/include/Exception_SetupOption.hpp"

using namespace std;

namespace OP2A{
namespace Setup{



template <typename T>
struct indent : public std::unary_function<T, void>
{
	indent(std::ostream& out) : os(out) {	}

	void operator()(T x)
	{
		os << "    " << x << "\n";
	}

	std::ostream& os;
};

const string SetupObject::NEST_SEPARATOR = ".";
const string SetupObject::NEST_SEPARATOR_ENVVAR = "_";


/*
 * Constructor
 */
//////////////////////////////////////////////////////////////////////////////
SetupObject::SetupObject(const string& name, const bool& helpOnError, const string& usageSummary)
	: NamedObject(name),
    m_options(),
    m_nested_objs(),
    m_setup(false),
    m_helpOnError(helpOnError),
    m_hadError(false),
    m_registered(false),
    m_setup_file(),
    m_env_vars(),
    m_nest(),
    m_usageSummary(usageSummary),
    parentSetup (OPNULL)
{
	op_assert(name != "");
}

/*
 * Virtual destruct
 */
SetupObject::~SetupObject()
{
	unregistFromSetupFacility();

	// removing object from parent vector
	if(parentSetup != OPNULL)
	{
		parentSetup->unregisterFromParent(this);
		parentSetup = OPNULL;
	}
}


/*
 * Public functions
 */
void SetupObject::unregistFromSetupFacility()
{
	if (m_registered)	SetupFacility::getInstance().unregistSetupObj(this);

	m_registered = false;
}




void SetupObject::setup (SetupArgs& args )
{
	m_setup = true;

	std::string basename = getNestName() + SetupObject::NEST_SEPARATOR;

	// prepend to each option the basename so they can be matched
	m_options.setHierarchyOptions (basename);

	// process the config arguments
	processOptions (args);

	// remove the base name form all options
	m_options.setHierarchyOptions ( std::string() );
}



void SetupObject::setupDynamic ( SetupArgs& args )
{
  std::string basename = getNestName() + SetupObject::NEST_SEPARATOR;

  // prepend to each option the basename so they can be matched
  m_options.setHierarchyOptions ( basename );

  // do dynamic configuration
  bool dynamic = true;
  m_options.processArgs ( args, dynamic );

  // remove the base name form all options
  m_options.setHierarchyOptions ( std::string() );
}


void SetupObject::setNest (const std::string& nest)
{
	m_nest = nest;
}

std::string SetupObject::getNest() const
{
	return m_nest;
}


std::string SetupObject::getNestName() const
{
	return m_nest.empty() ? getName() : m_nest +  SetupObject::NEST_SEPARATOR + getName();
}


void SetupObject::SetupNested (SetupObject* nestSetup, SetupArgs& args )
{
	op_assert(nestSetup != OPNULL);
	SetupNested(*nestSetup, args);
}


void SetupObject::SetupNested (SetupObject& nestSetup, SetupArgs& args )
{
  addNestedSetupObject(nestSetup);

  nestSetup.setNest (getNestName());      	// sets the name of the nesting object in this object
  nestSetup.setup(args);  					// do the configuration configure
}


void SetupObject::addNestedSetupObject(SetupObject& obj)
{
	m_nested_objs.push_back(&obj);

	std::sort(m_nested_objs.begin(), m_nested_objs.end(), std::less<SetupObject*>());
	std::vector<SetupObject*>::iterator last_obj = std::unique(m_nested_objs.begin(),m_nested_objs.end());
	m_nested_objs.erase(last_obj, m_nested_objs.end());

	obj.parentSetup = this;
}


void SetupObject::unsetup()
{
	m_setup = false;
}


std::string SetupObject::usageSummary() const
{
	return m_usageSummary;
}


std::string SetupObject::writeUsage() const
{
	std::ostringstream os;

	os << "Usage: " << getName() << " " << usageSummary() << "\n";
	os << "\n";
	os << "Options:\n";

	const OptionList& options = getOptionList();
	os << options.writeUsage();

	const std::vector<std::string>& env_vars = getEnvVars();
	if (!env_vars.empty())
	{
		os << "\nEnvironment Variables\n";
		for_each(env_vars.begin(), env_vars.end(), indent< std::string >(os));
	}

	std::string conffile ( getSetupFile());
	if (!conffile.empty())
	{
		os << "\nConfiguration File:" << conffile  << "\n";
	}

	return os.str();
}


Option * SetupObject::getOption(const std::string& name) const
{
	// gets the option
	Option * opt = m_options.findFirstByName(name);
	if (opt == OPNULL)	throw ExceptionSetupOption (FromHere(), "No option exists with name : " + name);

	return opt;
}



/*
 * Private Functions
 */
void SetupObject::unregisterFromParent(SetupObject * config)
{
	std::sort(m_nested_objs.begin(), m_nested_objs.end(), std::less<SetupObject*>());
	std::vector<SetupObject*>::iterator last_obj = std::unique(m_nested_objs.begin(),m_nested_objs.end());
	m_nested_objs.erase(last_obj, m_nested_objs.end());

	std::vector <SetupObject*>::iterator itr = find(m_nested_objs.begin(), m_nested_objs.end(), config);

	if ( itr != m_nested_objs.end() )
	{
		m_nested_objs.erase(itr);
	}
}






/*
 * Protect Functions
 */
void SetupObject::processOptions (SetupArgs& args)
{
	SetupArgs& parent_args = args;

	// default does not propagate file config out of the object
	SetupArgs file_args;
	processSetupFile(file_args);

	// default does not propagate envvars config out of the object
	SetupArgs env_args;
	processEnvironmentVariables (env_args);
	m_options.processArgs (env_args);

	m_options.doValidation();
	m_options.doComplete();
}

void SetupObject::processSetupFile (SetupArgs& args)
{
	if( !m_setup_file.empty() ) // only parse if filename was set
	{
		if (Common::fileExist(m_setup_file) != true)
		{
			throw Common::ExceptionFileSystem(FromHere(), "Could not open configuration file: " + m_setup_file);
		}

		SetupFileReader setup_parser;
		setup_parser.parse(m_setup_file, args);
  }

  m_options.processArgs ( args );
}

void SetupObject::processEnvironmentVariables (SetupArgs& args )
{
	std::pair<bool, std::string> ret;

	for (vector<std::string>::iterator it   = m_env_vars.begin(), stop = m_env_vars.end(); it != stop; ++it)
	{
		const std::string& envVar = *it;
		ret = readEnvironmentVariable(envVar);
		// The first member indicates if the second is valid or not
		if (ret.first) { args[envVar] = ret.second; }
	}

	m_options.processArgs ( args );
}

std::pair<bool, std::string>
SetupObject::readEnvironmentVariable(const std::string& envVar)
{
	std::string envVarNoDot(envVar);

	Common::StringOps::subst(NEST_SEPARATOR,NEST_SEPARATOR_ENVVAR, envVarNoDot);

	OPchar* ptrValue = getenv(envVarNoDot.c_str());

	if (ptrValue)
	{
		return std::pair<bool,std::string>(true, ptrValue);
	}
	else
	{
		return std::make_pair(false, std::string());
	}
}

void SetupObject::addEnvironmentVariable(const std::string& envVar)
{
	m_env_vars.push_back(envVar);
}

void SetupObject::registToSetupFacility()
{
	if (!m_registered)	SetupFacility::getInstance().registSetupObj(this);

	m_registered = true;
}


}
}
