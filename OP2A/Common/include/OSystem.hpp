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

#include "./Common/include/Common.hpp"

namespace OP2A
{
namespace Common{

class Common_API OSystem : public Common::NonCopyable<OSystem>
{


private:
	Common::ProcessInfo 	*m_process_info;
	Common::SignalHandler 	*m_sig_handler;
	Common::LibLoader 		*m_lib_loader;
};

} // end namespace Common
} // end namespace OP2A



#endif /* OSYSTEM_HPP_ */
