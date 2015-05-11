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

#include "Common/include/OP2A.hpp"




namespace OP2A{

namespace Common{

template <typename T>
class Array2D
{
public:
	Array2D():Size_(0), ElementSize_(0)
	{

	}


	Array2D(unsigned int I, unsigned int J):ElementSize_(J)
	{
		Resize(I);
	}


	// Getting information about the size
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


	// Resize data (I)
	void Resize(unsigned int I)
	{
		op_assert(ElementSize_ > 0);
		Data_.resize(I * ElementSize_);
		Size_	= I;
	}


	// Resize data (J)
	void SetElementSize (unsigned int J)
	{
		op_assert (J);
		const unsigned int OldSize = GetSize();
		ElementSize_ = J;
		Resize (OldSize);
	 }


	// Reserve data size
    void Reserve (unsigned int I)
    {
    	Data_.reserve (I*ElementSize_);
    }


    T & Get (unsigned int I, unsigned int J)
    {
    	op_assert (I < Size_);
    	op_assert (I < ElementSize_);
    	return Data_[I*ElementSize_+ J];
     }

    const T & Get (unsigned int I, unsigned int J) const
    {
    	op_assert (I < Size_);
        op_assert (I < ElementSize_);
        return Data_[I*ElementSize_+ J];
    }


    T & PushBack ()
    {
    	Resize (Size_+1);
    	return Get(Size_-1, 0);
    }

    void clear ()
    {
    	Resize (0);
    }



protected:
	unsigned int Size_;
	unsigned int ElementSize_;
	std::vector<T>	Data_;
};



}

}


#endif /* ARRAY2D_HPP_ */
