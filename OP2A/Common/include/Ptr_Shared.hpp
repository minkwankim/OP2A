/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Ptr_Shared.hpp
 * 			-  
 *  
 */
#ifndef PTR_SHARED_HPP_
#define PTR_SHARED_HPP_

#include "Common/include/Ptr_Safe.hpp"

namespace OP2A{
namespace Common{

/*
 * Class that defines a smart pointer that allows to share objects
 * (but ONLY objects deriving from OwnedObject) using reference counting.
 */

template<class T>
class Ptr_Shared
{
public:

	/// Constructor.
	Ptr_Shared ();
	explicit Ptr_Shared ( T* ptr );
	Ptr_Shared (const Ptr_Shared& other );
	Ptr_Shared (const Ptr_Shared<T>& other );

	/// Destructor.
	~Ptr_Shared();

	/// Release ownership
	void release();

	/// Reset the ptr
	void reset(T* other);

	/// Reset the ptr
	void reset(const Ptr_Shared& other);


	inline bool isNull() 	const { return (m_ptr == OPNULL); }
	inline bool isNotNull() const { return (m_ptr != OPNULL); }


	inline const Ptr_Shared& operator= (const Ptr_Shared& other);
	inline const Ptr_Shared& operator= (T* other);
	inline bool operator== (const Ptr_Shared& other);


	T* getPtr() const { return m_ptr; }

	T* operator->() const
	{
		op_assert(m_ptr != OPNULL);
		return m_ptr;
	}

	T& operator*() const
	{
		op_assert(m_ptr != OPNULL);
		return *m_ptr;
	}

private:
	T* m_ptr;	// Raw pointer to the object
};


/*
 * Definition of functions
 */

// Constructors
template<class T>
Ptr_Shared<T>::Ptr_Shared() : m_ptr(0) {}

template<class T>
Ptr_Shared<T>::Ptr_Shared(T* ptr)
{
	m_ptr = ptr;
	if (m_ptr != OPNULL)	m_ptr->addOwner();
}

template<class T>
Ptr_Shared<T>::Ptr_Shared(const Ptr_Shared& other) : m_ptr(other.m_ptr)
{
  if (m_ptr != OPNULL) {
    m_ptr->addOwner();
  }
}

template<class T>
Ptr_Shared<T>::Ptr_Shared(const Ptr_Shared<T>& other) : m_ptr(other.m_ptr)
{
  if (m_ptr != OPNULL)
    m_ptr->addOwner();
}


// Destructor
template<class T>
Ptr_Shared<T>::~Ptr_Shared()
{
	if(m_ptr != OPNULL)
	{
		m_ptr->removeOwner();
		if(m_ptr->hasNoOwner())	deletePtr(m_ptr);
	}
}



template<class T>
void Ptr_Shared<T>::release()
{
	if (m_ptr != OPNULL)
	{
		m_ptr->removeOwner();
		if(m_ptr->hasNoOwner())	deletePtr(m_ptr);
	}

	m_ptr = OPNULL;
}


template<class T>
void Ptr_Shared<T>::reset(T* other)
{
	if (m_ptr != OPNULL)
	{
		m_ptr->removeOwner();
		if(m_ptr->hasNoOwner())	deletePtr(m_ptr);
	}

	m_ptr = other;
	if (m_ptr != OPNULL)	m_ptr->addOwner();
}

template<class T>
void Ptr_Shared<T>::reset(const Ptr_Shared& other)
{
	if (m_ptr != OPNULL)
	{
		m_ptr->removeOwner();
		if(m_ptr->hasNoOwner())	deletePtr(m_ptr);
	}

	m_ptr = other.m_ptr;
	if (m_ptr != OPNULL)	m_ptr->addOwner();
}



// Operators
template<class T>
inline const Ptr_Shared<T>& Ptr_Shared<T>::operator= (const Ptr_Shared& other)
{
	reset(other);
	return *this;
}

template<class T>
inline const Ptr_Shared<T>& Ptr_Shared<T>::operator= (T* other)
{
	reset(other);
	return *this;
}

template<class T>
inline bool Ptr_Shared<T>::operator== (const Ptr_Shared& other)
{
  return (m_ptr == other.m_ptr);
}


}
}



#endif /* PTR_SHARED_HPP_ */
