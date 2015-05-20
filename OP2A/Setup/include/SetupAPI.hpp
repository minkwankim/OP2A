/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 18, 2015
 *      			Author: Minkwan Kim
 *
 * SetupAPI.hpp
 * 			-  
 *  
 */
#ifndef SETUPAPI_HPP_
#define SETUPAPI_HPP_

#include "Common/include/ExportAPI.hpp"


/*
 * Define the macro for Setup_API
 *
 * @author	Minkwan Kim
 * @version 1.0	18/5/2015
 */

#ifdef Setup_EXPORT
#	define Setup_API	OP_EXPORT_API
#else
#	define Setup_API		OP_IMPORT_API
#endif


namespace OP2A{
namespace Setup{

}
}

#endif /* SETUPAPI_HPP_ */
