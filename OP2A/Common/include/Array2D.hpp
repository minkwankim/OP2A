/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 11, 2015
 *      			Author: Minkwan Kim
 *
 * Array2D.hpp
 * 			-  
 *  
 */
#ifndef ARRAY2D_HPP_
#define ARRAY2D_HPP_

#include "Common/include/Standard_headers.hpp"
#include "Common/include/OPAssert.hpp"



namespace OP2A{

namespace Common{

/*
 * Class for 2D array (Matrix, I x J)
 * 	- [Size_]: size of Coloumn (I)
 * 	- [ElementSize_]: Size of row (J)
 * 	- [Data_]: Data of matrix
 *
 * 	Initially written by:	Minkwan Kim
 * 	Last modified on	:	May/14/2015
 * 	Last modified by	:	Minkwan Kim
 */

template <typename T>
class Array2D
{
protected:
	unsigned int Size_;
	unsigned int ElementSize_;
	std::vector<T>	Data_;

public:

	// Constructors
	Array2D():Size_(0), ElementSize_(0)
	{

	}

	Array2D(unsigned int _I, unsigned int _J):ElementSize_(_J)
	{
		Resize(_I);
	}




	/*
	 * Class Functions
	 */
	// CF-01: Getting information about the size
    unsigned int GetSize () const
    {
    	return Size_;
    };

    unsigned int size () const
    {
    	return GetSize();
    }

    unsigned int GetElementSize () const
    {
    	return ElementSize_;
    }


	// CF-02: Resize data (I)
	void Resize(unsigned int I)
	{
		op_assert(ElementSize_ > 0);
		Data_.resize(I * ElementSize_);
		Size_	= I;
	}

	// CF-03: Resize data (J)
	void SetElementSize (unsigned int J)
	{
		op_assert (J);
		const unsigned int OldSize = GetSize();
		ElementSize_ = J;
		Resize (OldSize);
	 }


	// CF-04: Reserve data size
    void Reserve (unsigned int I)
    {
    	Data_.reserve (I*ElementSize_);
    }


    // CF-05: Get data
    T & Get (unsigned int I, unsigned int J)
    {
    	op_assert (I < Size_);
    	op_assert (J < ElementSize_);
    	return Data_[I*ElementSize_+ J];
     }

    const T & Get (unsigned int I, unsigned int J) const
    {
    	op_assert (I < Size_);
        op_assert (I < ElementSize_);
        return Data_[I*ElementSize_+ J];
    }


	// CF-04: Push back data
    T & PushBack ()
    {
    	Resize (Size_+1);
    	return Get(Size_-1, 0);
    }

    // CF-05: Clear data
    void clear ()
    {
    	Resize (0);
    }
};



}

}


#endif /* ARRAY2D_HPP_ */
