/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Ptr_Safe.hpp
 * 			-  
 *  
 */
#ifndef PTR_SAFE_HPP_
#define PTR_SAFE_HPP_


#include "Common/include/Exception_FailedCast.hpp"
#include "Common/include/Ptr_Alloc.hpp"
//#include "Common/include/DemangledTypeID.hpp"

namespace OP2A{
namespace Common{

template < class TYPE > class Ptr_Shared;

template <class TYPE>
class Ptr_Safe
{
	friend class Ptr_Shared<TYPE>;

public:
	Ptr_Safe(): m_ptr(OPNULL) {};
	Ptr_Safe(TYPE *ptr){  m_ptr = ptr; };
	Ptr_Safe(const Ptr_Safe& other){ m_ptr = other.m_ptr;};

	~Ptr_Safe(){};

	/// Reset the ptr
	void reset(TYPE* other){ m_ptr = other;};
	void reset(const Ptr_Safe& other){ m_ptr = other.m_ptr;};

	/// Dynamic casts the m_ptr to the TYPE
	template < typename DTYPE > Ptr_Safe<DTYPE> d_castTo() const;


	/// Check if it is CFNULL
	inline bool isNull() const;

	/// Check if it is Not CFNULL
	inline bool isNotNull() const;

	/// Overloading of "=" with SafePtr
	inline const Ptr_Safe& operator= (const Ptr_Safe& other);
	inline const Ptr_Safe& operator= (TYPE* other);

	/// Overloading of "=="
	inline bool operator== (const Ptr_Safe& other) const;

	/// Overloading of "!="
	inline bool operator!= (const Ptr_Safe& other) const;

	/// Overloading of "<"
	inline bool operator< (const Ptr_Safe& other) const;

	/// Overloading of ">"
	inline bool operator> (const Ptr_Safe& other) const;

	/// Overloading of "*"
	TYPE& operator*() const {return *m_ptr;}

	/// Overloading of "->"
	TYPE* operator->() const {return m_ptr;}

private:
	TYPE *m_ptr;


public:
	template < typename Ret >
    class mem_fun_t
    {
    public:
		typedef TYPE	argument_type;
		typedef Ret		result_type;

	public:
		/// constructor
		explicit mem_fun_t( Ret ( TYPE::*pf )() ) : m_f(pf) {}
		Ret operator()(const Ptr_Safe<TYPE>& p) const  { return (p.m_ptr->*m_f)(); }

	private: // data
		Ret (TYPE::*m_f)();
    };
};


template < typename TYPE , typename Ret >
inline typename Common::Ptr_Safe<TYPE>::template mem_fun_t<Ret>
safeptr_mem_fun( Ret (TYPE::*f)() )
{
	return typename Common::Ptr_Safe<TYPE>::template mem_fun_t<Ret> (f);
}




template< typename TYPE>
inline bool Ptr_Safe<TYPE>::isNull() 	const { return (m_ptr == OPNULL); }

template< typename TYPE>
inline bool Ptr_Safe<TYPE>::isNotNull() const { return !isNull(); }

template< typename TYPE>
inline const Ptr_Safe<TYPE>& Ptr_Safe<TYPE>::operator= (const Ptr_Safe& other)
{
  reset(other);
  return *this;
}

template< typename TYPE>
inline const Ptr_Safe<TYPE>& Ptr_Safe<TYPE>::operator= (TYPE* other)
{
  reset(other);
  return *this;
}

template< typename TYPE>
inline bool Ptr_Safe<TYPE>::operator== (const Ptr_Safe& other) const
{
  return (m_ptr == other.m_ptr);
}

template< typename TYPE>
inline bool Ptr_Safe<TYPE>::operator!= (const Ptr_Safe& other) const
{
  return (m_ptr != other.m_ptr);
}

template< typename TYPE>
inline bool Ptr_Safe<TYPE>::operator< (const Ptr_Safe& other) const
{
  return (m_ptr < other.m_ptr);
}

template< typename TYPE>
inline bool Ptr_Safe<TYPE>::operator> (const Ptr_Safe& other) const
{
  return (m_ptr > other.m_ptr);
}


template< typename TYPE>
template< typename DTYPE>
Ptr_Safe<DTYPE> Ptr_Safe<TYPE>::d_castTo() const
{
  DTYPE* rPtr = dynamic_cast<DTYPE*>(m_ptr);
  if (rPtr == OPNULL)
  {
    std::string msg ("SafePtr failed dynamic cast from ");
    //msg += DEMANGLED_TYPEID(TYPE);
    msg += " to ";
   // msg += DEMANGLED_TYPEID(DTYPE);
    throw Common::ExceptionFailedCast(FromHere(), msg);
  }

  return Ptr_Safe<DTYPE>(rPtr);
}


}
}





#endif /* PTR_SAFE_HPP_ */
