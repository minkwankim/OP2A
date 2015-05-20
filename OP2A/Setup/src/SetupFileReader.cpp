/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 18, 2015
 *      			Author: Minkwan Kim
 *
 * SetupFileReader.hpp
 * 			-  
 *  
 */

#include <fstream>

#include "Common/include/OP2A.hpp"
#include "Common/include/StringOps.hpp"

#include "Common/include/Exception_FileSystem.hpp"
#include "Setup/include/SetupFileReader.hpp"
#include "Setup/include/SetupFileReader_ASCII.hpp"

namespace OP2A{
namespace Setup{


void SetupFileReader::parse(const std::string& name, SetupArgs& args)
{
	// check the file type
	std::ifstream	setupfile(name.c_str());

	if(!setupfile) throw Common::ExceptionFileSystem (FromHere(), "Could not open file: " + name);

	std::string line;
	getline(setupfile, line);

	bool is_xml = Common::StringOps::startsWith(line, "<?xml");
	setupfile.close();


	SetupFileReader_ASCII setupFile;
	setupFile.parse(name, args);
}


}
}
