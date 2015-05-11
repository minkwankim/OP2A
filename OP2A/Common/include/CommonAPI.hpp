/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 11, 2015
 *      			Author: Minkwan Kim
 *
 * CommonAPI.hpp
 * 			-  
 *  
 */
#ifndef COMMONAPI_HPP_
#define COMMONAPI_HPP_

#include "Common/include/ExportAPI.hpp"

/*
 * Define the macro Common_API
 * 	- note build system defines Common_EXPORTS when compiling Common files
 */

#ifdef Common_EXPORTS
	#define Common_API	OP_EXPORT_API
#else
	#define Common_API	OP_IMPORT_API
#endif



#endif /* COMMONAPI_HPP_ */
