/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 15, 2015
 *      			Author: Minkwan Kim
 *
 * Time_HourMinSec.cpp
 * 			-  
 *  
 */

#include <sstream>
#include "Common/include/Time_HourMinSec.hpp"



namespace OP2A{
namespace Common{


const OPdouble	TimeData::secPerHour	= 1.0/3600.0;
const OPdouble	TimeData::secPerMin		= 1.0/60.0;
const OPdouble	TimeData::usecPerSec	= 1.0e-6;



TimeData::TimeData(): m_h(0.0), m_m(0.0), m_s(0.0)
{

}

TimeData::TimeData(const OPdouble &sec)
{
	set(sec);
}


TimeData::~TimeData()
{

}



void TimeData::set(double t_in_sec)
{
	m_h = floor(t_in_sec * secPerHour);
	m_m = floor((t_in_sec - m_h * 3600.) * secPerMin);
	m_s = t_in_sec - m_h * 3600. - m_m * 60.;
}


std::string TimeData::str() const
{
	std::ostringstream	oss;
	if (m_h > 0)	oss << m_h << " h ";
	if (m_m > 0)	oss << m_m << "min ";
	oss << m_s << "sec";

	return oss.str();
}


}
}
