/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 15, 2015
 *      			Author: Minkwan Kim
 *
 * Time_HourMinSec.hpp
 * 			-  
 *  
 */
#ifndef TIME_HOURMINSEC_HPP_
#define TIME_HOURMINSEC_HPP_

#include "Common/include/Common.hpp"
#include "Common/include/OP2A.hpp"
#include "Common/include/CommonAPI.hpp"



namespace OP2A{
namespace Common{

/*
 * Data for time
 *
 * @author	Minkwan
 * @version 1.0 15/5/2015
 */
class Common_API	TimeData{
public:
	TimeData();
	TimeData(const OPdouble &sec);

	~TimeData();

	std::string str() const;
	void set(const OPdouble t_in_sec);	// set the value correspondint to the give time in secconds

	static const OPdouble secPerHour;
	static const OPdouble secPerMin;
	static const OPdouble usecPerSec;

private:
	OPdouble m_h;
	OPdouble m_m;
	OPdouble m_s;
};


}
}



#endif /* TIME_HOURMINSEC_HPP_ */
