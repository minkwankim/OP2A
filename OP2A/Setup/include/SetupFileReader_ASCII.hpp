/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 18, 2015
 *      			Author: Minkwan Kim
 *
 * SetupFileReader_nested.hpp
 * 			-  
 *  
 */
#ifndef SETUPFILEREADER_NESTED_HPP_
#define SETUPFILEREADER_NESTED_HPP_


#include "Common/include/OP2A.hpp"

#include "Setup/include/SetupAPI.hpp"
#include "Setup/include/SetupArgs.hpp"


namespace OP2A{
namespace Setup{

/*
 * Parses a setup file in ASCII format
 * @author	Minkwan Kim
 * @version	1.0	18/5/2015
 */

class Setup_API SetupFileReader_ASCII {
public:
	/// Reads the given file into a set of name/value pairs, both as std::strings.
	explicit SetupFileReader_ASCII();

	/// Virtual destructor.
	virtual ~SetupFileReader_ASCII();

	/// Copys all of the name/value pairs into the supplied map of
	virtual void parse (const std::string& filename,	SetupArgs& args);

protected:
	std::string expandFileName(const std::string& fname) const;


private:
	// Indicates to continue in nextline.
	static const OPchar CONTINUATOR;

	// Separates tags and values.
	static const std::string SEPARATOR;

	// Initiates comments. Rest of the line is ignored.
	static const std::string COMMENTOR;

	// Initiates meta-comments
	static const std::string METACOMMENTOR;

	/// Initiates block comments. Starting from that point
	static const std::string BLOCKCOMMENTORBEGIN;
	static const std::string BLOCKCOMMENTOREND;
}; // Class NestedConfigFileReader

}
}

#endif /* SETUPFILEREADER_NESTED_HPP_ */
