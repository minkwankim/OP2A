/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * WorkerStatus.hpp
 * 			-  
 *  
 */
#ifndef WORKERSTATUS_HPP_
#define WORKERSTATUS_HPP_

#include "Common/include/Common.hpp"


namespace OP2A{
namespace Common{


class Common_API WorkerStatus
{
public:
	enum Type
	{
		INVALID		= -1,
		STARTING	= 0,
		CONFIGURING	= 1,
		EXITING		= 2,
		RUNNING		= 3,
		NOT_RUNNING	= 4,
		WAITING		= 5,
		IDLE		= 6
	};

	//typedef Common::EnumT<WorkerStatus>	ConvertBase;

	//struct

};


} // end namespace Common
} // end namepsace OP2A



#endif /* WORKERSTATUS_HPP_ */
