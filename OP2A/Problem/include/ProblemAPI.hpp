/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 21, 2015
 *      			Author: Minkwan Kim
 *
 * ProblemAPI.hpp
 * 			-  
 *  
 */
#ifndef PROBLEMAPI_HPP_
#define PROBLEMAPI_HPP_

#include "Common/include/ExportAPI.hpp"


#ifdef Problem_EXPORTS
#	define Problem_API			OP_EXPORT_API
#   define Problem_TEMPLATE
#else
#   define Problem_API			OP_IMPORT_API
#   define Problem_TEMPLATE		OP_TEMPLATE_EXTERN
#endif


/// Macro defining the factories in Environment
#ifdef  HAVE_CXX_EXPLICIT_TEMPLATES
#	define Problem_Factory(__fac__) namespace OP2A { namespace Environment { OP_TEMPLATE_EXTERN template class Problem_API OP2A::Problem::Factory<OP2A::Problem::__fac__>; } }
#else
#   define Problem_Factory(__fac__)
#endif



#endif /* PROBLEMAPI_HPP_ */
