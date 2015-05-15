/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 11, 2015
 *      			Author: Minkwan Kim
 *
 * OPAssert.hpp
 * 			-  
 *  
 */
#ifndef OPASSERT_HPP_
#define OPASSERT_HPP_



#include <cassert>
#include "Common/include/NonCopyable.hpp"
#include "Common/include/CommonAPI.hpp"


namespace OP2A{

/// Manager of behavior of assertions
class Common_API AssertionManager : public Common::NonCopyable <AssertionManager>
{
public:
	  bool DoAssertions;		// Flag to dynamically turn off assertions
	  bool AssertionDumps;		// Flag for Assertion Dumps backtraces
	  bool AssertionThrows;		// Flag for Assettion Throws exceptions

	  /// Constructor
	  AssertionManager();

	  /*
	   * Class functions
	   */

	  // [CF-01] Gets the instance of the manager
	  static	AssertionManager&	getInstance ();

	  // [CF-02] If AssertionManager is not handling assertions and those are passed to the standard assert function
	  // Controlled by the build option ENABLE_STDASSERT
	  static bool notHandlingAssertions ()
	  {
	  	  #ifdef ENABLE_STDASSERT
		  	  return true;
	  	  #else
		  	  return false;
	  	  #endif
	  }

	  // [CF-03] Forward declaration of the function that implements the always present assert
	  static void do_assert (bool condition, const char * cond_str, const char * file, int line, const char * func, const char * desc = 0 );
};






#ifndef ENABLE_STDASSERT
#define op_always_assert(a) \
     { if (!(a)) { ::OP2A::AssertionManager::do_assert((a), #a, __FILE__, __LINE__, __FUNCTION__); } }

#define op_always_assert_desc(msg, a) \
     { if (!(a)) { ::OP2A::AssertionManager::do_assert((a), #a, __FILE__, __LINE__, __FUNCTION__, msg); } }


/// Assertions are off if compiled with DNDEBUG
#ifndef NDEBUG
	#define op_assert(a)       		op_always_assert((a))
	#define op_assert_desc(m,a) 	op_always_assert_desc(m,(a))
#else
	#define op_assert(a)
	#define op_assert_desc(m,a)
#endif
#else
//ENABLE_STDASSERT
#ifndef NDEBUG
	#define op_assert(a)              	 assert(a)
	#define op_assert_desc(m,a)          assert(a)
	#define op_always_assert(a)          assert(a)
	#define op_always_assert_desc(msg,a) assert(a)
#else
	#define op_assert(a)
	#define op_assert_desc(m,a)
	#define op_always_assert(a)
	#define op_always_assert_desc(msg,a)
#endif

#endif // CF_ENABLE_STDASSERT

//////////////////////////////////////////////////////////////////////////////

}


#endif /* OPASSERT_HPP_ */
