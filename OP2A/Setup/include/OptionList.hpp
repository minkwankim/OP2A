/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 26, 2015
 *      			Author: Minkwan Kim
 *
 * OptionList.hpp
 * 			-  
 *  
 */
#ifndef OPTIONLIST_HPP_
#define OPTIONLIST_HPP_

#include "Setup/include/Option.hpp"
#include "Setup/include/OptionT.hpp"
#include "Setup/include/SetupArgs.hpp"
#include "Setup/include/Exception_BadMatch.hpp"


namespace OP2A{
namespace Setup{


/*
 * Class to process and to store the list of options
 */
class Setup_API OptionList {

public: // functions

	/*
	 * Constructor
	 */
	OptionList();						// Options are added after construction.
	OptionList(const OptionList& obj);	// Copy constructor

	/// Destructor
	~OptionList();


	/// Assignement operator
	const OptionList& operator=(const OptionList& obj);

	/// Processes the supplied argument list.
	/// The arguments that are not "consumed" can be retrieved via
	/// <code>getUnprocessedArgs</code>.
	/// @param args Arguments to be processed
	/// @param dynamic indicates if the arguments come from a dynamic source, for example interactive input
	/// @return a pair with the processed arguments and the unprocessed args
	std::pair<SetupArgs, SetupArgs> processArgs (SetupArgs& args, bool dynamic = false );

	/// Process command line options
	void processCommandLine(int argc, char** argv);

	/// Writes the usage of each option to the given output stream.
	/// The format is that as for --help:<p>
	/// <pre>
	/// --opt     Option description (default = default)
	/// </pre>
	/// Note that there is one tab between the option and the description.
	/// Also note that there is no "Option:" line, since these lists can be chained together.
	std::string writeUsage() const;

	/// Writes all the options, including those in the configuration file.
	std::string writeConfigUsage() const;

	/// Returns an XML with the option descriptions
	/// @returns string with xml
	std::string getOptionsXML() const;

	/// Adds the optio to the this list.
	/// Note that ownership of the option object is taken by this object.
	/// @param opt the pointer to the option to add. Will be deleted when OptionList is destroyed.
	/// @returns the option that was added.
	Option* addOptionPtr(Option* opt);

	/// Makes the first option take on the value from the second option,
	///  if it is not set explicitly by the user
	/// @param from name of option to take value from
	/// @param to name of option to set value
	void linkOptions(const std::string& from, const std::string& to);

	/// Nests all the options in the correct nest label.
	/// @param std::string with the nested label.
	void setHierarchyOptions(const std::string& nest);

	/// Prepend all the options in the supplied string if they do not start with it
	/// @param std::string with the prepending name
	void prependToName(const std::string& pre);

	///  Returns the options for the given name.
	/// @param name full name of the Option
	std::vector<Option*> findByName(const std::string& name) const;

	/// Returns the first option for the given name.
	/// @returns the Option pointer if found or CFNULL if not found
	/// @see findByName
	Option* findFirstByName(const std::string& name) const;


	/// Checks for existence of an option with the same name in this optionvset.
	/// @throws DuplicateNameException when matches are found.
	void checkNameDoesNotExist(const std::string& name) const;




	/// Gets the option list debug info
	/// @returns string with info
	std::string debugInfo() const;



	/// Sets if command line arguments passed can have semantical mistakes
	void setStrictArgs(const bool& strict) { m_strictArgs = strict; }

	/// Gets if command line arguments passed can have semantical mistakes
	bool getStrictArgs() const { return m_strictArgs; }


	/// Fires the onComplete method on all options
	void doComplete();

	/// Fires the validate method on all options
	void doValidation();



	/// Adds the option, associated with the given label name
	/// @param name label of the option
	/// @param description description of the option
	template < typename Type >
	Option* addConfigOption(const std::string& name,
	const std::string& description)
	{
		OptionT<Type> * ptr = new OptionT<Type>(name, description);
		return addOptionPtr(ptr);
	}

	/// Adds the option, associated with the given label name
	/// Marks the option using a Trait class with type MarkOption
	/// @param name label of the option
	/// @param description description of the option
	template < typename Type, typename MarkOption >
	Option* addConfigOption(const std::string& name,
	const std::string& description)
	{
		OptionT<Type> * ptr = new OptionT<Type>(name, description);
		MarkOption::mark(ptr);
		return addOptionPtr(ptr);
	}

	/// Adds the option, associated with the given label name
	/// @param name label of the option
	/// @param description description of the option
	template < typename OptionType >
	Option* addConfigOptionByType(const std::string& name,
	const std::string& description)
	{
		OptionType * ptr = new OptionType(name, description);
		return addOptionPtr(ptr);
	}

	/// Adds the option, associated with the given label name
	/// Marks the option using a Trait class with type MarkOption
	/// @param name label of the option
	/// @param description description of the option
	template < typename OptionType, typename MarkOption >
	Option* addConfigOptionByType(const std::string& name,
	const std::string& description)
	{
		OptionType * ptr = new OptionType(name, description);
		MarkOption::mark(ptr);
		return addOptionPtr(ptr);
	}












private:
	/// Removes all Option's from the list
	void deleteOptions();

	/// Returns the names of all options that match the label.
	/// @param label string against which to compare the options label
	std::vector<std::string> findMatches(const std::string& label) const;

	/// Complains about the number of matches then throws an exception.
	/// @throws DuplicateNameException when multiple matches are found
	void complainAboutMatches(const std::vector<std::string>& matches, const std::string& tag) const;



private: // member data
	/// The option handlers registered for this set.
	std::vector<Option*> m_vec_options;

	/// Storage of the options that are dynamic
	std::vector<Option*> m_dyn_options;

	/// command line arguments passed cannot have semantical mistakes
	bool m_strictArgs;

}; // class OptionList





}
}



#endif /* OPTIONLIST_HPP_ */
