/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * SignalHandler.hpp
 * 			-  
 *  
 */
#ifndef SIGNALHANDLER_HPP_
#define SIGNALHANDLER_HPP_

#include "Common/include/NonCopyable.hpp"
#include "Common/include/Object_Owned.hpp"


namespace OP2A{
namespace Common{

class Common_API SignalHandler : public Common::NonCopyable<SignalHandler>
{
public:

	/// signal function type
	typedef void (*sighandler_t)(int);

	SignalHandler();
	virtual ~SignalHandler();

	/// Regists the signal handlers that will be handled by this class
	virtual void registSignalHandlers() = 0;
	static std::string getClassName() { return "SignalHandler"; }

protected: // methods

private: // data
};

}
}

#endif /* SIGNALHANDLER_HPP_ */
