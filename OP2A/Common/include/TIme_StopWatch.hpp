/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 15, 2015
 *      			Author: Minkwan Kim
 *
 * TIme_StopWatch.hpp
 * 			-  
 *  
 */
#ifndef TIME_STOPWATCH_HPP_
#define TIME_STOPWATCH_HPP_

#include "Common/include/Time_HourMinSec.hpp"
#include "Common/include/Time_Info.hpp"


namespace OP2A{
namespace Common{


/*
 * Measure elapsed seconds
 *
 * @author	Minkwan Kim
 * @version 1.0 15/5/2015
 */
template <typename TIMEDATA = Common::WallTime>
class Time_StopWatch
{
private:
	bool		m_running;
	OPdouble	m_total;
	TIMEDATA	impl;

private:
	void initStartTime();
	void takeStopTime();
	void accumulateTime();

public:
	Time_StopWatch():m_running(false), m_total(0.0) { };
	~Time_StopWatch() { };



	// Starting timming from 0.0
	// @post	m_running = true
	// @post	m_total = 0.0
	void start();

	/// Restart timing from 0.00. Same as calling start()
	/// @post m_running == true
	/// @post m_total == 0
	void restart();

	/// Reset the timer
	/// Clears the elapsed time.
	/// @post m_running == false
	/// @post m_total == 0
	void reset();

	/// Resumes counting. Doesn't clear the elapsed time.
	/// No effect if isRunning()
	/// @post m_running == true
	/// @post m_total >= 0
	void resume();

	/// Stop timing. Doesn't clear the elapsed time.
	/// No effect if isNotRunning()
	/// @post m_running == false
	/// @post m_total >= 0
	void stop();

	/// Read current time.
	/// @return current time in seconds
	OPdouble read() const;

	/// Converst the current time to CFdouble
	/// @return current time in seconds
	operator OPdouble () const { return read(); }

	/// Read current time.
	/// @return current time in Hour, Minutes and Seconds
	TimeData readTimeHMS();

	/// Checks if the Stopwatch is m_running
	/// @return TRUE if it is running
	bool isRunning() const;

	/// Checks if the Stopwatch is not running
	/// @return TRUE if it is not running
	bool isNotRunning() const;
};





template <typename TIMEDATA>
void Time_StopWatch<TIMEDATA>::start()
{
	reset();
	resume();
}


template <typename TIMEDATA>
void Time_StopWatch<TIMEDATA>::restart()
{
	start();
}


template <typename TIMEDATA>
void Time_StopWatch<TIMEDATA>::reset()
{
	stop();
	m_total = 0.;
}


template <typename TIMEDATA>
void Time_StopWatch<TIMEDATA>::resume()
{
	if (isNotRunning())
	{
		initStartTime();
		m_running = true;
	}
}


template <typename TIMEDATA>
void Time_StopWatch<TIMEDATA>::stop()
{
	if (isRunning())
	{
		takeStopTime();
		accumulateTime();
		m_running = false;
	}
}


template <typename TIMEDATA>
OPdouble Time_StopWatch<TIMEDATA>::read() const
{
	if (isNotRunning())	return m_total;

	return impl.getDeltaT() + m_total;
}


template <typename TIMEDATA>
TimeData Time_StopWatch<TIMEDATA>::readTimeHMS()
{
	return TimeData(read());
}


template <typename TIMEDATA>
bool Time_StopWatch<TIMEDATA>::isRunning() const
{
	return m_running;
}


template <typename TIMEDATA>
bool Time_StopWatch<TIMEDATA>::isNotRunning() const
{
	return !isRunning();
}


// Private functions
template <typename TIMEDATA>
void Time_StopWatch<TIMEDATA>::initStartTime()
{
	  impl.initStartTime();
}

template <typename TIMEDATA>
void Time_StopWatch<TIMEDATA>::takeStopTime()
{
	  impl.takeStopTime();
}

template <typename TIMEDATA>
void Time_StopWatch<TIMEDATA>::accumulateTime()
{
	  impl.accumulateTime();
}


}
}



#endif /* TIME_STOPWATCH_HPP_ */
