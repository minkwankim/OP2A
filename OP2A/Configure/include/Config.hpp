/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 11, 2015
 *      			Author: Minkwan Kim

 *	[This is adapted from CoolFluid of VKI] 			-
 *
 */

#ifndef CONFIGURE_CONFIG_HPP
#define CONFIGURE_CONFIG_HPP


#include "Common/include/ExportAPI.hpp"


/// Define the macro Config_API
/// @note build system defines Config_EXPORTS when compiling Config files
#ifdef CONFIG_EXPORTS
#	define CONFIG_API OP_EXPORT_API
#else
#	define CONFIG_API OP_IMPORT_API
#endif




namespace OP2A {
namespace Config {

}
}

#endif // COOLFluiD_Config_hh
