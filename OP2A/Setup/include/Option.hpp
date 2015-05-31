/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 26, 2015
 *      			Author: Minkwan Kim
 *
 * Option.hpp
 * 			-  
 *  
 */
#ifndef OPTION_HPP_
#define OPTION_HPP_

#include "Common/include/OP2A.hpp"
#include "Common/include/Ptr_Alloc.hpp"

#include "Setup/include/SetupAPI.hpp"


namespace OP2A{
namespace Setup{




class OptionList;
class OptionValidation;

class Setup_API Option
{
public:

	/*
	 * Constructors
	 */
	Option(const std::string& label, const std::string& description); 	/// Constructor using a label and a description
	Option(const Option& rhs);											/// Copy constructor


	/*
	 * Virtual destructor
	 */
	virtual ~Option();


	/// Returns the label of this option.
	virtual std::string getLabel() const;

	/// Sets the label of this option.
	virtual void setLabel(const std::string& label);

	/// Returns the nested label of this option.
	virtual std::string getNestedLabel() const;

	/// Sets the hierarchy owner of this option.
	virtual void setHierarchyOwner(const std::string& owner);

	/// Returns the hierarchy owner of this option.
	virtual std::string getHierarchyOwner() const;

	/// Returns the description.
	virtual std::string description() const;

	/// Links an external owned value to the Option
	/// @param value pointer to the parameter to link
	virtual void linkToValue(void* value, const std::string& typestr) = 0;

	/// Sets the default value
	/// @param defvalue pointer to the value to set
	virtual void setDefaultValue(void* defvalue, const std::string& typestr) = 0;

	/// Sets this option value from a std::string
	/// @pre label has already been checked to ensure that it matches the value's tag
	/// @param value string with value to be converted to actual type
	virtual void setValue(const std::string& value) = 0;

	/// Sets all linked options to the given value.
	virtual void setLinkedOptions(const std::string& value);

	/// Returns the default value, as a std::string
	virtual std::string getDefaultValueAsString() const = 0;

	/// Returns the value, as a std::string
	virtual std::string getValueAsString() const = 0;

	/// Gets the option description in a --help style format.
	virtual std::string helpDescription() const;

	/// Called when the option is successfully processed.
	/// By default it does nothing and returns true.
	/// Subclasses should reimplement to add specific actions.
	/// @return true if successful
	virtual void onProcessed();

	/// Called when the option is updated via an interactive change
	/// By default it does nothing and returns true.
	/// Subclasses should reimplement to add specific actions.
	virtual void onUpdate();

	/// Called after all options have been successfully processed.
	/// By default it does nothing and returns true.
	/// Subclasses should reimplement to add specific actions.
	virtual void onComplete(const OptionList& opts);

	/// Returns whether the name matches the nested label
	virtual bool matchNestedLabel(const std::string& label) const;

	/// Adds the given option as being tied to this one.
	/// If this option is set, so too will the other option.
	virtual void addLinkedOption(Option* const opt);

	/// Test a string to see if it mathces the common delimiter
	static bool matchesDelimiter(const std::string& arg);

	/// Checks the argument has a dilimiter
	/// @throw BadMatchException if doens't have a dilimiter
	void checkArgument(const std::string& arg);

	/// Validate this option by using all pointers to OptionValidation
	/// @throw OptionValidationException when one of the validations fails
	virtual void validate ();

	/// Adds the given validation to this option
	/// @param val pointer to the validation
	/// @post passed validation will be owned by this option
	virtual void addValidation(OptionValidation* val);

	/// Makes this option interactive
	virtual void makeOptionDynamic();

	/// Marks this option as basic
	virtual void makeOptionBasic();

	/// Checks if the option is interactive
	bool isDynamic() const
	{
		return m_dynamic;
	}

	/// Checks if the option is basic
	bool isBasic() const
	{
		return m_basic;
	}


	/// Returns the low-level details about this object
	/// @return string with details
	virtual std::string debugInfo() const = 0;

	/// Returns XML description of this option
	/// @return string with details
	virtual std::string getXML () const = 0;

	/// Sets values from command line parameters
	virtual OPuint setCommandLineValues(const std::string& arg, OPuint position, int argc, char** argv) = 0;

	/// Creates a new Option* with the same values as this
	virtual Option* clone() const = 0;


public: // data
	static const std::string DELIMITER;


protected: // methods
	/// Returns the options that are linked to this one
	/// @return vector with pointer to options
	std::vector<Option*>& getLinkedOptions();


private: // methods
	/// deletes all pointers to validations
	void clearValidations();


private: // data
	/// Options that inherit from this one.
	std::vector<Option*> m_linked_opts;

	/// list of option validations
	std::vector<OptionValidation*> m_validations;

	/// The label of this option
	std::string m_label;

	/// The description of this option
	std::string m_description;

	/// The nested label for this option
	std::string m_owner_hname;

	/// Is this option dynamic
	bool m_dynamic;

	/// Is this option basic
	bool m_basic;

}; // class Option




}
}



#endif /* OPTION_HPP_ */
