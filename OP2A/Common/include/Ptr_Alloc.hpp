/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Ptr_Alloc.hpp
 * 			-  
 *  
 */
#ifndef PTR_ALLOC_HPP_
#define PTR_ALLOC_HPP_

#include "Common/include/OPAssert.hpp"


namespace OP2A{


#define OPNULL	0

template <class TYPE>
void deletePtr(TYPE*& ptr)
{
	if (ptr != OPNULL)
	{
		delete ptr;
		ptr	= OPNULL;
	}

	op_assert(ptr == OPNULL);
}

template <class TYPE>
void deletePtrAtrray(TYPE *&ptr)
{
	if (ptr != OPNULL)
	{
		delete [] ptr;
		ptr = OPNULL;
	}

	op_assert(ptr == OPNULL);
}

}



#endif /* PTR_ALLOC_HPP_ */
