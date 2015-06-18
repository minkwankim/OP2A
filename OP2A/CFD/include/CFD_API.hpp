/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 16, 2015
 *      			Author: Minkwan Kim
 *
 * CFD_API.hpp
 * 			-  
 *  
 */
#ifndef CFD_API_HPP_
#define CFD_API_HPP_


#include "Common/include/ExportAPI.hpp"


/*
 * Define the macro for CFD_API
 *
 * @author	Minkwan Kim
 * @version 1.0	14/6/2015
 */

#ifdef CFD_EXPORT
#	define CFD_API	OP_EXPORT_API
#else
#	define CFD_API	OP_IMPORT_API
#endif

namespace OP2A{
namespace CFD{

}
}



#endif /* CFD_API_HPP_ */
