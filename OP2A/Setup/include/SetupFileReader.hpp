/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 18, 2015
 *      			Author: Minkwan Kim
 *
 * SEtupFileReader.hpp
 * 			-  
 *  
 */
#ifndef SETUPFILEREADER_HPP_
#define SETUPFILEREADER_HPP_

#include "Common/include/OP2A.hpp"

#include "Setup/include/SetupAPI.hpp"
#include "Setup/include/SetupArgs.hpp"


namespace OP2A{
namespace Setup{


/*
 * Class for setup file reading
 * @author	Minkwan Kim
 * @version	1.0	20/05/2015
 */

class Setup_API SetupFileReader
{
public:
	SetupFileReader()	{	};
	~SetupFileReader()	{	};

	void parse(const std::string& name, SetupArgs& args);
};







}
}




#endif /* SETUPFILEREADER_HPP_ */
