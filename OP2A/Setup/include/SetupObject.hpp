/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 21, 2015
 *      			Author: Minkwan Kim
 *
 * SetuoObject.hpp
 * 			-  
 *  
 */
#ifndef SETUOOBJECT_HPP_
#define SETUOOBJECT_HPP_


#include "Common/include/NamedObject.hpp"

#include "Setup/include/SetupAPI.hpp"
#include "Setup/include/OptionList.hpp"




namespace OP2A{
namespace Setup{

class Setup_API	SetupObject : public Common::NamedObject
{
public: // functions

	/// Constructor
	SetupObject(const std::string& name, const bool& helpOnError = false, const std::string& usageSummary = std::string());

	/// Virtual destructor
	virtual  ~SetupObject();


	/*
	 * Public functions
	 */
	/// Removes this object registration from the SetupFacility
	void unregistFromSetupFacility();

	/// Setup the options for this object.
	/// @param args is the SetupArgs with the arguments to be parsed.
	virtual void setup ( SetupArgs& args );

	/// Setup the dynamic options for this object.
	/// To be extended by derived classes.
	/// @param args is the SetupArgs with the arguments to be parsed.
	virtual void setupDynamic ( SetupArgs& args );

	/// Sets this object _nest member.
	/// @param nest the name of the parent object that nests this one
	void setNest(const std::string& nest);

	/// Gets this object _nest member.
	/// @return string with nest label
	std::string getNest() const;

	// Gets this object full nested name.
	/// @return _nest + "." + _objName
	std::string getNestName() const;

	/// Setup the nested objects of this object with the arguments passed onto the object.
	virtual void SetupNested(SetupObject* nestSetup, SetupArgs& args );
	virtual void SetupNested(SetupObject& nestSetup, SetupArgs& args );

	/// Adds a nested SetupObject
	void addNestedSetupObject(SetupObject& obj);

	/// Unsetups the object
	virtual void unsetup();


	/// Check if the object is Setupured
	virtual bool isSetup() const
	{
		return m_setup;
	}


	/// Provides a usage summary for all options of the object
	virtual std::string usageSummary() const;


	/// Returns the usage of this SetupObject
	std::string writeUsage() const;

	/// Gets an option by its name
	Option * getOption(const std::string& name) const;


	/// Sets the Setup options by calling the defineSetupOptions
	template <typename CLASS>
	void addSetupOptionsTo(const CLASS* ptr)
	{
		CLASS::defineSetupOptions(m_options);
		registToSetupFacility();
	}

	/// Gets an option by its name, but on its derived type
	template < typename OPTIONTYPE >
	OPTIONTYPE * getOptionT(const std::string& name) const
	{
		Option * opt = getOption(name);
		OPTIONTYPE * dptr = dynamic_cast<OPTIONTYPE*>(opt);
		if (dptr == OPNULL)
		{
			std::string msg ("Option failed dynamic cast to ");
			msg += DEMANGLED_TYPEID(OPTIONTYPE);
			throw Common::ExceptionFailedCast(FromHere(),msg);
		}
		return dptr;
	}


	/// Associates the specified parameter to an option identified by its name.
	template < typename TYPE ,  typename DEFTYPE  >
	void setParameter(const std::string& name, TYPE* const value, DEFTYPE defvalue)
	{
		Option * opt = getOption(name);

		// links the value to the option
		opt->linkToValue(value, DEMANGLED_TYPEID(TYPE));

		// sets the default value
		opt->setDefaultValue(&defvalue, DEMANGLED_TYPEID(DEFTYPE));
	}


	OptionList getOptionList() const
	{
		return m_options;
	}


	OptionList& getOptionList()
	{
		return m_options;
	}


	std::vector<std::string> getEnvVars() const
	{
		return m_env_vars;
	}

	std::string getSetupFile() const
	{
		return m_setup_file;
	}


	virtual void setSetupFile (const std::string file)
	{
		m_setup_file = file;
	}


	/// Returns whether there was an error during processing.
	virtual bool hadError() const
	{
		return m_hadError;
	}


protected:
	/// Processes the argument list from the Setup arguments
	/// Derived classes may override to change order in which options are parsed
	/// Default is:
	///   1) Setup file via processSetupFile()
	///   2) parent arguments via processParentArgs()
	///   3) environmental variables via processEnvironmentVariables()
	/// Last Setup has preference
	virtual void processOptions (SetupArgs& args);

	/// Parses and processes the Setupu file
	virtual void processSetupFile (SetupArgs& args);

	/// Processes the options defined by the environmental variables
	virtual void processEnvironmentVariables (SetupArgs& args);


	/// Reads the environment variable, splits it by whitespace and processes it as a command line.
	virtual std::pair<bool,std::string> readEnvironmentVariable(const std::string& envVar);

	/// Adds an environment variable to be parsed.
	virtual void addEnvironmentVariable(const std::string& envVar);

	/// Regists this object to the SetupFacility
	void registToSetupFacility();


private: // functions
	/// Sets the SetupObject given by parameter as the parent of this SetupObject
	//void setParentSetup(SetupObject * parentSetup);

	/// Unregisters (erases from the vector) the given child SetupObject object.
	void unregisterFromParent(SetupObject * Setup);



	/*
	 * Member data
	 */

public:
	/// Separator for nested Setups
	static const std::string NEST_SEPARATOR;

	/// Separator for nested Setup but for the enviromental variables
	static const std::string NEST_SEPARATOR_ENVVAR;

private:
	OptionList m_options;						/// the list of options of this object
	std::vector< SetupObject* > m_nested_objs;	/// the list of nested Setup objects

	bool m_setup;								/// Setup flag to indicate that the object is configurated.
	bool m_helpOnError;							/// Whether to show help when processing fails.
	bool m_hadError;							/// Whether there was an error when processing the arguments.
	bool m_registered;							/// Whether the object has been registered to the SetupFacility

	std::string m_setup_file;					/// filename for setup

	std::vector<std::string> m_env_vars;		/// The environment variables.
	std::string m_nest;							/// The nested names of this object.
	std::string m_usageSummary;					/// Summary of usage.

	SetupObject * parentSetup;					/// Pointer to the parent SetupObject
};






}
}


#endif /* SETUOOBJECT_HPP_ */
