/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * ProcessInfo.hpp
 * 			-  
 *  
 */
#ifndef PROCESSINFO_HPP_
#define PROCESSINFO_HPP_

#include "Common/include/Common.hpp"
#include "Common/include/NonCopyable.hpp"
#include "Common/include/Compatibility.hpp"
#include "Common/include/OP2A.hpp"

namespace OP2A{
namespace Common{

/*
 * Class to representation of the current information on the memory usage.
 *
 * 	Initially written by:	Minkwan Kim
 * 	Last modified on	:	May/14/2015
 * 	Last modified by	:	Minkwan Kim
 */


class Common_API ProcessInfo :public Common::NonCopyable<ProcessInfo>
{
public:
	// Constructro and dectructor
	ProcessInfo();			// Constructor
	virtual ~ProcessInfo();	// Destructor

	/*
	 * Class Functions
	 */
	// CF-01: For Name of Flatform
	virtual std::string getFlatformName()	const = 0;

	// CF-02: Dump BAcktrace
	virtual std::string getBackTrace()		const = 0;

	// CF-03: Process ID
	virtual OPuint 		getPID () 			const = 0;

	// CF-04: The memory usage in bytes
	virtual OPdouble memoryUsageBytes()	const = 0;

	// CF-05: String with the memory usage
	std::string memoryUsage()			const;

	// CF-06: Class name
	static std::string getClassName()
	{
		return "ProcessInfo";
	}
};

} //end namespace Common
} // end namepsace OP2A






#endif /* PROCESSINFO_HPP_ */
