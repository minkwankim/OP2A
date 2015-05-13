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

class Common_API ProcessInfo :public Common::NonCopyable<ProcessInfo>
{
public:

	ProcessInfo();	// Constructor
	virtual ~ProcessInfo();	// Destructor

	virtual std::string getFlatformName()	const = 0;
	virtual std::string getBackTrace()		const = 0;
	virtual OPuint 		getPID () 			const = 0;

	virtual OPdouble memoryUsageBytes()	const = 0;
	std::string memoryUsage()			const;
	static std::string getClassName()
	{
		return "ProcessInfo";
	}
};

} //end namespace Common
} // end namepsace OP2A






#endif /* PROCESSINFO_HPP_ */
