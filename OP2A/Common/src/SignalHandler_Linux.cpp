/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * SignalHandler_Linux.cpp
 * 			-  
 *  
 */


#include <cstdio>     // for printf()
#include <cstdlib>    // for free() and abort()
#include <csignal>    // POSIX signal(), SIGFPE and SIGSEGV
#include <fenv.h>     // floating Common access
#include <sstream>    // streamstring

#include "Common/include/Common.hpp"
#include "Common/include/Exception_FloatingPoint.hpp"
#include "Common/include/ProcessInfo_Linux.hpp"
#include "Common/include/SignalHandler_Linux.hpp"



namespace OP2A{
namespace Common{



SignalHandlerLinux::SignalHandlerLinux()
{
}

SignalHandlerLinux::~SignalHandlerLinux()
{
}



void SignalHandlerLinux::registSignalHandlers()
{
	// register handler functions for the signals
	signal(SIGFPE,    (sighandler_t) SignalHandlerLinux::handleSIGFPE);
	signal(SIGSEGV,   (sighandler_t) SignalHandlerLinux::handleSIGSEGV);

	// enable the exceptions that will raise the SIGFPE signal
	feenableexcept ( FE_DIVBYZERO );
	feenableexcept ( FE_INVALID   );
	feenableexcept ( FE_OVERFLOW  );
	feenableexcept ( FE_UNDERFLOW );
}


int SignalHandlerLinux::handleSIGFPE (int signal)
{
	printf("\nreceived signal SIGFPE [%d] - 'Floating Point Exception'\n",signal);
	static std::string dump = ProcessInfoLinux::dumpBacktrace();
	printf( "%s\n", dump.c_str() );
	throw Common::ExceptionFloatingPoint(FromHere(), "Some floating point operation has given an invalid result");
}


int SignalHandlerLinux::handleSIGSEGV(int signal)
{
	printf("\nreceived signal SIGSEGV [%d] - 'Segmentation violation'\n",signal);
	static std::string dump = ProcessInfoLinux::dumpBacktrace();
	printf( "%s\n", dump.c_str() );
	abort();
}


}
}
