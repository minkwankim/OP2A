/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * SignalHandler_Linux.hpp
 * 			-  
 *  
 */
#ifndef SIGNALHANDLER_LINUX_HPP_
#define SIGNALHANDLER_LINUX_HPP_

#include "Common/include/SignalHandler.hpp"


namespace OP2A{
namespace Common{

/// This class handles of signals from the Linux operating system
class Common_API SignalHandlerLinux : public SignalHandler {

public:

  SignalHandlerLinux();
  virtual ~SignalHandlerLinux();

  virtual void registSignalHandlers();

protected: // methods
  static int handleSIGFPE(int signal);
  static int handleSIGSEGV(int signal);
};


}
}


#endif /* SIGNALHANDLER_LINUX_HPP_ */
