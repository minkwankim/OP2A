/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * ProcessInfo_Linux.cpp
 * 			-  
 *  
 */
#ifndef PROCESSINFO_LINUX_HPP_
#define PROCESSINFO_LINUX_HPP_



#include "Common/include/ProcessInfo.hpp"


namespace OP2A{
namespace Common{

class Common_API ProcessInfoLinux :public ProcessInfo
{
public:

	ProcessInfoLinux();	// Constructor
	virtual ~ProcessInfoLinux();	// Destructor

	virtual std::string getFlatformName()	const	{ return "Linux";};
	static 	std::string dumpBacktrace();
	virtual std::string getBackTrace()		const;

	virtual OPuint 		getPID () 			const;
	virtual OPdouble 	memoryUsageBytes()	const;
	std::string memoryUsage()				const;
};

} //end namespace Common
} // end namepsace OP2A




#endif /* PROCESSINFO_LINUX_CPP_ */
