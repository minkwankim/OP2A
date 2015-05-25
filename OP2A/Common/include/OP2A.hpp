/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 11, 2015
 *      			Author: Minkwan Kim
 *
 * OPPA.hpp
 * 			-  
 *  
 */
#ifndef OPPA_HPP_
#define OPPA_HPP_





#include "Common/include/Compatibility.hpp"
#include "Common/include/Standard_headers.hpp"
#include "Common/include/Common.hpp"
#include "Common/include/OPAssert.hpp"


#include "Common/include/Version.hpp"







/* Definition of OP2A namespace
 * 			by Minkwan Kim
 */
namespace  OP2A {

	/// Definition of the basic types for possible portability conflicts
	typedef float              OPfloat;
	typedef double             OPdouble;
	typedef long double        OPldouble;


#ifdef HAVE_LONG
	typedef long long int      		OPint;
	typedef unsigned long long int	OPuint;
#else
	typedef int                		OPint;
	typedef unsigned int       		OPuint;
#endif

	typedef char						OPchar;
	typedef std::string					OPstring;
	//typedef std::complex				OPcomplex;		// Need to fix later

	/// Enumeration of the dimensions
	enum Dim					{DIM_0D, DIM_1D, DIM_2D, DIM_3D};

	/// Enumeration of the coordinates indexes
	enum CoordXYZ				{XX, YY, ZZ};

	/// Enumeration of the reference coordinates indexes
	enum CoordRefXiEtaZeta    {KSI, ETA, ZTA};

	/// Enumeration of the device types
	enum DeviceType 			{CPU=0, GPU=1, MIC=2};


	/// class to be used to define a default type
	class NOTYPE {};


	///////////////////////////////////////
	/// Definition of the default precision
	///////////////////////////////////////
	#ifdef PRECISION_LONG_DOUBLE
		typedef OPldouble OPreal;
	#else
		#ifdef PRECISION_DOUBLE
			typedef OPdouble OPreal;
		#else
			#ifdef PRECISION_SINGLE
				typedef OPfloat OPreal;
			#endif
		#endif
	#endif

	// if nothing defined, use double
	#if !defined PRECISION_DOUBLE && !defined PRECISION_SINGLE && !defined PRECISION_LONG_DOUBLE
		typedef OPdouble OPreal;
	#endif
}
// namespace OP2A


#endif /* OPPA_HPP_ */
