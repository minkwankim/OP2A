/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * OSystem.cpp
 * 			-  
 *  
 */


#include "Common/include/OSystem.hpp"
#include "Common/include/ProcessInfo.hpp"
#include "Common/include/SignalHandler.hpp"
#include "Common/include/Ptr_Safe.hpp"
#include "Common/include/StringOps.hpp"

#ifdef OS_LINUX
	#include "Common/include/ProcessInfo_Linux.hpp"
	#include "Common/include/SignalHandler_Linux.hpp"
#endif

#ifdef OS_MACOSX
	#include "Common/include/ProcessInfo_OSX.hpp"
	#include "Common/include/SignalHandler_OSX.hpp"
#endif

#ifdef OS_WINDOWS
	#include "Common/include/ProcessInfo_WIN.hpp"
	#include "Common/include/SignalHandler_WIN.hpp"
	#include "Common/include/Win32LibLoader.hpp"
#endif


namespace OP2A{
namespace Common{


OSystem::OSystem():
		m_process_info(OPNULL),
		m_sig_handler(OPNULL),
		m_lib_loader(OPNULL)
{

#ifdef OS_LINUX
    if (m_process_info == OPNULL )	m_process_info = new ProcessInfoLinux();
    if ( m_sig_handler == OPNULL )  m_sig_handler = new SignalHandlerLinux();
#else
	#ifdef OS_MACOSX
    	if ( m_process_info == OPNULL ) m_process_info = new ProcessInfoMacOSX();
    	if ( m_sig_handler == OPNULL )  m_sig_handler = new SignalHandlerMacOSX();
	#else
		#ifdef OS_WINDOWS
    		if ( m_process_info == OPNULL ) m_process_info = new ProcessInfoWin32();
    		if ( m_sig_handler == OPNULL )  m_sig_handler = new SignalHandlerWin32();
    		if ( m_lib_loader == OPNULL )   m_lib_loader = new Win32LibLoader();
		#else
  	  	  #error	"Unkown operating system: not Windows, MacOSX or Linux"
		#endif
	#endif
#endif
}


OSystem::~OSystem()
{
	deletePtr(m_process_info);
	deletePtr(m_sig_handler);
	deletePtr(m_lib_loader);
}



OSystem& OSystem::getInstance()
{
	static OSystem osystem;
	return osystem;
}

Ptr_Safe<ProcessInfo> OSystem::getProcessInfo()
{
	return m_process_info;
}


Ptr_Safe<SignalHandler> OSystem::getSignalHandler()
{
	return m_sig_handler;
}

Ptr_Safe<LibLoader> OSystem::getLibLoader()
{
	return m_lib_loader;
}


void OSystem::executeCommand(const std::string& call)
{
	int return_value = system ( call.c_str() );
	if ( return_value == -1)
	{
		std::string msg;
		msg += "Command \'";
		msg += call;
		msg += "\' return error code";
		throw ExceptionOSystem ( FromHere(), msg );
	}
}


void OSystem::sleep (const OPuint& seconds)
{
	std::string callSleep = "sleep " + Common::StringOps::to_str(seconds);
	executeCommand (callSleep);
}



}
}
