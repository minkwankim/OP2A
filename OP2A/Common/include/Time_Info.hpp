/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 15, 2015
 *      			Author: Minkwan Kim
 *
 * TimeInfo.hpp
 * 			-  
 *  
 */
#ifndef TIMEINFO_HPP_
#define TIMEINFO_HPP_


#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "Common/include/Standard_headers.hpp"
#include "Common/include/OP2A.hpp"
#include "Common/include/CommonAPI.hpp"


//////////////////////////////////////////////////////////////////////////////
namespace OP2A {
namespace Common{

/*
 * Class for Time usages check
 *
 * 	@author	Minkwan Kim
 * 	@version 1.0	15/5/2015
 */
class Common_API TimeinfoRusage
{
private:
	OPdouble m_StartTime;	// Starting time
	OPdouble m_StopTime;	// Stopping time

private:
	OPdouble seconds() const;	// Get time difference

public:
	TimeinfoRusage();

	void initStartTime();
	void takeStopTime();

	// @return Elapsed time in seconds
	OPdouble getDeltaT() const;

	void accumulateTime(OPdouble& accTime);
};

class Common_API TimeinfoCclock
{
private:
	OPdouble m_StartTime;	// Starting time
	OPdouble m_StopTime;	// Stopping time

private:
	OPdouble seconds() const;	// Get time difference

public:
	TimeinfoCclock();

	void initStartTime();
	void takeStopTime();

	// @return Elapsed time in seconds
	OPdouble getDeltaT() const;

	void accumulateTime(OPdouble& accTime);
};



typedef TimeinfoRusage CPUTime;
typedef TimeinfoCclock WallTime;
}
}
#endif /* TIMEINFO_HPP_ */
