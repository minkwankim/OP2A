/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 11, 2015
 *      			Author: Minkwan Kim
 *
 * Compatability.hpp
 * 			-  
 *  
 */
#ifndef COMPATABILITY_HPP_
#define COMPATABILITY_HPP_



namespace OP2A {

/// Macro necessary to cope with operator<< in combination with templates.
/// Compilers that have to have this set: cxx ( SGI )
/// Compilers that should not have this set: gcc, CC, icc (8.0)
#ifdef CXX_NEEDS_FRIEND_TMPL_DECL
	#define LTGT <>
#else
	#define LTGT
#endif

/// Macro necessary for compilers, that do not define the __FUNCTION__ variable
#ifndef HAVE_FUNCTION_DEF
	#define __FUNCTION__ ""
#endif

/// Macro necessary for compiling functions to be used on device and/or host
#ifndef HAVE_CUDA
	#define HOST_DEVICE
#else
	#define HOST_DEVICE __host__ __device__
#endif

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif /* COMPATABILITY_HPP_ */
