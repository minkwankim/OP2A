/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 15, 2015
 *      			Author: Minkwan Kim
 *
 * TimeInfo.cpp
 * 			-  
 *  
 */



#include <time.h>

#include "Common/include/Time_Info.hpp"
#include "Common/include/Time_HourMinSec.hpp"

//////////////////////////////////////////////////////////////////////////////
namespace OP2A {
namespace Common{


/*
 * Class for Time usages check
 *
 * 	@author	Minkwan Kim
 * 	@version 1.0	15/5/2015
 */
// Constructor
TimeinfoRusage::TimeinfoRusage() : m_StartTime(0.0), m_StopTime(0.0) { }


OPdouble TimeinfoRusage::seconds() const
{
	rusage usg;
	getrusage(RUSAGE_SELF, &usg);
	timeval usertime = usg.ru_utime;
	return (usertime.tv_sec + (usertime.tv_usec * TimeData::usecPerSec));
}


void TimeinfoRusage::initStartTime()
{
	m_StartTime	= seconds();
}

void TimeinfoRusage::takeStopTime()
{
	m_StopTime	= seconds();
}

void TimeinfoRusage::accumulateTime(OPdouble& accTime)
{
	accTime	+= (m_StopTime - m_StartTime);
}

OPdouble TimeinfoRusage::getDeltaT() const
{
	return (seconds() - m_StartTime);
}








/*
 * Class for Time usages check using cclock
 *
 * 	@author	Minkwan Kim
 * 	@version 1.0	15/5/2015
 */
// Constructor
TimeinfoCclock::TimeinfoCclock() : m_StartTime(0.0), m_StopTime(0.0) { }


OPdouble TimeinfoCclock::seconds() const
{
	const OPdouble	secs_per_tick = 1.0 / CLOCKS_PER_SEC;
	return ( static_cast<OPdouble>(clock()) ) * secs_per_tick;
}


void TimeinfoCclock::initStartTime()
{
	m_StartTime	= seconds();
}

void TimeinfoCclock::takeStopTime()
{
	m_StopTime	= seconds();
}

void TimeinfoCclock::accumulateTime(OPdouble& accTime)
{
	accTime	+= (m_StopTime - m_StartTime);
}

OPdouble TimeinfoCclock::getDeltaT() const
{
	return (seconds() - m_StartTime);
}








}
}
