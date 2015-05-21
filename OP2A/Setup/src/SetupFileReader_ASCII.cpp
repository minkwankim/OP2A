/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 18, 2015
 *      			Author: Minkwan Kim
 *
 * SetupFileReader_ASCII.cpp
 * 			-  
 *  
 */


#include <fstream>

#include "Common/include/CodeLocation.hpp"
#include "Common/include/StringOps.hpp"
#include "Common/include/Exception_FileSystem.hpp"

#include "Setup/include/SetupArgs.hpp"
#include "Setup/include/SetupFileReader.hpp"
#include "Setup/include/SetupFileReader_ASCII.hpp"
#include "Setup/include/Exception_BadMatch.hpp"



using namespace  std;

namespace OP2A{
namespace Setup
{

const std::string	SetupFileReader_ASCII::SEPARATOR   = "=";
const OPchar		SetupFileReader_ASCII::CONTINUATOR = '\\';
const std::string	SetupFileReader_ASCII::COMMENTOR   = "//";
const std::string	SetupFileReader_ASCII::METACOMMENTOR   = "///";
const std::string	SetupFileReader_ASCII::BLOCKCOMMENTORBEGIN   = "/*";
const std::string	SetupFileReader_ASCII::BLOCKCOMMENTOREND   = "*/";


SetupFileReader_ASCII::SetupFileReader_ASCII()
{
}

SetupFileReader_ASCII::~SetupFileReader_ASCII()
{
}




// Get expanded file
std::string SetupFileReader_ASCII::expandFileName(const std::string& fname) const
{
	string	exp	= fname;
	if (fname[0] == '~')
	{
		char* home	= getenv("HOME");

		if (home)
		{
			exp	= string(home)	+ fname.substr(1);
		}
	}

	return exp;
}







void SetupFileReader_ASCII::parse (const std::string& filename_in, SetupArgs& args)
{
	if(!filename_in.empty())
	{
		std::string filename = expandFileName(filename_in);	// Get filename

		std::ifstream iss(filename.c_str());				// Open file
		if(!iss) throw Common::ExceptionFileSystem (FromHere(),"Could not open file: " + filename);



		// Initialize comment block indicator
		bool commentedBlock = false;

		while (iss)
		{
			// Get line from the file
			std::string line;
			getline(iss, line);

			// 1. Removing commented blocks
			if(commentedBlock == true)
			{
				std::string::size_type blockendcmtpos = line.find(BLOCKCOMMENTOREND);
				if (blockendcmtpos != std::string::npos)
				{
					// strip to the end of the line
					commentedBlock = false;
					line.erase(0,blockendcmtpos + BLOCKCOMMENTOREND.size());
				}
				else
				{
					line.erase();
				}
			}


			// 2. Check for meta-comment with include directive
			std::string::size_type metacmtpos = line.find(METACOMMENTOR);
			if ( metacmtpos != std::string::npos )
			{
				std::vector<std::string> wds = Common::StringOps::getWords(line);

				if (wds.size() > 1 ) // more than metacomment is present
				{
					// check for include directive
					if ( wds[1] == "IncludeCase" )
					{
						if ( wds.size() > 2 )
						{
							SetupFileReader other_file_reader;
							other_file_reader.parse(wds[2], args);
						}
						else
						{
							//CFLogWarn ( "IncludeCase directive has no file defined\n" );
						}
					}
				}
			}


			// 3. Check for start of a commented block
			std::string::size_type blockbegincmtpos = line.find(BLOCKCOMMENTORBEGIN);

			if (blockbegincmtpos != std::string::npos)
			{
				line.erase(blockbegincmtpos);		// strip to the end of the line
				commentedBlock = true;
			}


			/// 4. Removing comments
			std::string::size_type cmtpos = line.find(COMMENTOR);
			if (cmtpos != std::string::npos)
			{
				line.erase(cmtpos);
			}


			// Parse the read line
			Common::StringOps::trimRear(line);	// strip trailing whitespaces

			if(!line.empty())
			{
				// if the last character is a backslash, keep concatenating the lines
				while (iss && line[line.length() - 1] == CONTINUATOR)
				{
					line = std::string(line, 0, line.length() - 1);

					std::string nextline;
					getline(iss, nextline);

					Common::StringOps::trimRear(nextline);
					line += nextline;
				}

				OPint	shift = 0;
				size_t	eqpos = line.find(SEPARATOR);

				if (eqpos != std::string::npos)
				{
					std::string preSeparator(line, eqpos-1, 1);

					if(preSeparator == "+")
					{
						shift = 1;
					}

					std::string label(line, 0, eqpos-shift);
					std::string value(line, eqpos + 1);

					Common::StringOps::trim(label); // trim front and back around the label
					Common::StringOps::trim(value); // they might have leading and tailing whitespace



					if((preSeparator == "+") && (args[label] != ""))
					{
						args[label] = args[label] + " " + value;
					}
					else
					{
						if(args[label] != "")
						{
							throw Setup::ExceptionBadMatch(FromHere(),"Overriding previously specified option for: " + label);
						}
						args[label] = value;
					}
				}
			} // not empty

		} // while loop
	} // filename is empty
}








}
}


