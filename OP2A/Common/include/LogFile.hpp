/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 20, 2015
 *      			Author: Minkwan Kim
 *
 * LogFile.hpp
 * 			-  
 *  
 */
#ifndef LOGFILE_HPP_
#define LOGFILE_HPP_

#include	"Common/include/OP2A.hpp"


namespace OP2A{

enum	LogLevel
{
	EMERG	= 0,
	FATAL	= 0,
	ALERT	= 100,
	CRIT	= 200,
	ERROR	= 300,
	WARN	= 400,
	NOTICE	= 500,
	INFO	= 600,
	DEBUG	= 700,
	NOTSET	= 800
};


class Common_API	Logger : public Common::NonCopyable<Logger>
{
public:
	static	Logger& getInstance();


private:
	Logger();
	~Logger();
};


#define OPout    Logger::getInstance().getMainLogger().noticeStream()
#define OPerr    Logger::getInstance().getMainLogger().errorStream()
#define OPlog    Logger::getInstance().getMainLogger().infoStream()
#define OPtrace  Logger::getInstance().getTraceLogger().debugStream()

#define LogInfo(x)   OPlog(INFO,x)
#define LogNotice(x) OPlog(NOTICE,x)
#define LogWarn(x)   OPlog(WARN,x)


}




#endif /* LOGFILE_HPP_ */
