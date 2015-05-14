/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 14, 2015
 *      			Author: Minkwan Kim
 *
 * MemoryAllocator.hpp
 * 			-  
 *  
 */
#ifndef MEMORYALLOCATOR_HPP_
#define MEMORYALLOCATOR_HPP_

#include "Common/include/Exception_MemoryAllocator.hpp"


namespace OP2A{
namespace Common{

/*
 * Class for Memory Allocation
 * 	- This class uses virtual functions, because the overhead is neglectable here.
 *
 * 	Initially written by:	Minkwan Kim
 * 	Last modified on	:	May/14/2015
 * 	Last modified by	:	Minkwan Kim
 */
class Common_API MemoryAllocator
{
public:
	typedef void * MA_Ptr;
	typedef size_t MA_Size;

	virtual ~MemoryAllocator () { }



	/*
	 * Class Functions
	 */
	// CF-01: Return current size (bytes)
	virtual MA_Size GetSize () const = 0;

	// CF-02: Granularity
	virtual MA_Size GetGranularity () const = 0;

	// CF-03: Resizes memory and Returns new size
	virtual MA_Size Resize (MA_Size NewSize) = 0;

	// CF-04: Is state valid (= memory allocated?)
	virtual bool IsValid () const = 0;

	// CF-05: Return pointer to memory
	virtual MA_Ptr GetPtr () const = 0;

	// CF-06: Is resize zero copy?
	virtual bool IsZeroCopy () const = 0;


	operator MA_Ptr () const
	{
		return GetPtr();
	};

	operator bool () const
	{
		return IsValid ();
	};
};


}
}


#endif /* MEMORYALLOCATOR_HPP_ */
