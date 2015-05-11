/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 11, 2015
 *      			Author: Minkwan Kim
 *
 * ExportAPI.hpp
 * 			-  
 *  
 */
#ifndef EXPORTAPI_HPP_
#define EXPORTAPI_HPP_



// Define or not the extern keyword for templates
#ifdef HAVE_EXTERN_TEMPLATES
	#define OP_TEMPLATE_EXTERN extern
#else
	#define OP_TEMPLATE_EXTERN
#endif


// Visual Studio :
//    all symbols are local (invisible) by default
#if defined (_MSC_VER)
	#define OP_LOCAL_API
	#define OP_IMPORT_API    __declspec(dllimport)
	#define OP_EXPORT_API    __declspec(dllexport)
#endif

// GNU Compiler :  all symbols are global (visible) by default
#if defined (__GNUC__) // && defined (__unix__)
	#define OP_LOCAL_API    __attribute__ ((visibility("hidden")))
	#define OP_EXPORT_API   __attribute__ ((visibility("default")))
	#define OP_IMPORT_API   __attribute__ ((visibility("default")))
#endif


// For unrecognized compilers
#ifndef OP_LOCAL_API
	#if defined (ACCEPT_ANY_COMPILER)
		#define OP_LOCAL_API
		#define OP_EXPORT_API
		#define OP_IMPORT_API
	#else
		#error "Unrecognized compiler and/or platform."
	#endif
#endif



#endif /* EXPORTAPI_HPP_ */
