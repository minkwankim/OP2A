/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * OSystem.hpp
 * 			-  
 *  
 */
#ifndef OSYSTEM_HPP_
#define OSYSTEM_HPP_

#include "Common/include/Common.hpp"
#include "Common/include/NonCopyable.hpp"
#include "Common/include/Exception.hpp"
#include "Common/include/Exception_OSystem.hpp"
#include "Common/include/Ptr_Safe.hpp"
#include "Common/include/ProcessInfo.hpp"
#include "Common/include/SignalHandler.hpp"
#include "Common/include/LibLoader.hpp"


namespace OP2A
{
namespace Common{

class Common_API OSystem : public Common::NonCopyable<OSystem>
{


private:
	Common::ProcessInfo 	*m_process_info;
	Common::SignalHandler 	*m_sig_handler;
	Common::LibLoader 		*m_lib_loader;

	OSystem();
	~OSystem();

public:
	static OSystem& getInstance();

	Common::Ptr_Safe<Common::ProcessInfo> 	getProcessInfo();
	Common::Ptr_Safe<Common::SignalHandler> getSignalHandler();
	Common::Ptr_Safe<Common::LibLoader> 	getLibLoader();

	void executeCommand (const std::string& call);
	void sleep (const OPuint& seconds = 1);
};

} // end namespace Common
} // end namespace OP2A



#endif /* OSYSTEM_HPP_ */
