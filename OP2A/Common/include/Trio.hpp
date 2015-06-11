/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 27, 2015
 *      			Author: Minkwan Kim
 *
 * Trio.hpp
 * 			-  
 *  
 */
#ifndef TRIO_HPP_
#define TRIO_HPP_



#include "Common/include/OP2A.hpp"

namespace OP2A{
namespace Common{




template <typename T1, typename T2, typename T3>
class Trio
{
public:
	typedef T1	first_type;
	typedef T2	second_type;
	typedef T3	third_type;

	T1 first;
	T2 second;
	T3 third;

public:
	Trio()	: first(T1()), second(T2()), third(T3())
	{

	}

	Trio(const T1& a, const T2& b, const T3& c) : first(a), second(b), third(c)
	{

	}

	template <typename U1, typename U2, typename U3>
	Trio(const Trio<U1, U2, U3>& p) : first(p.first), second(p.second), third(p.third)
	{
	}


	friend std::istream& operator>> (std::istream& in, Trio<T1, T2, T3>& q)
	{
		in >> q.first  >> q.second >> q.third;
		return in;
	}

	friend std::ostream& operator<< (std::ostream& out, Trio<T1, T2, T3>& q)
	{
		out << q.first << " " << q.second << " " << q.third;
		return out;
	}
};


template <class T1, class T2, class T3>
inline bool operator==(const Trio<T1, T2, T3>& x, const Trio<T1, T2, T3>& y)
{
	return x.first == y.first && x.second == y.second && x.third == y.third;
}

template <class T1, class T2, class T3>
inline bool operator!=(const Trio<T1, T2, T3>& x, const Trio<T1, T2, T3>& y)
{
	return !(x == y);
}

template <class T1, class T2, class T3>
inline Trio<T1, T2, T3> make_Trio(const T1& x, const T2& y, const T3& z)
{
	return Trio<T1, T2, T3>(x, y, z);
}


template <typename T1, typename T2, typename T3> std::istream& operator>> (std::istream& in, Trio<T1, T2, T3>&);
template <typename T1, typename T2, typename T3> std::ostream& operator<< (std::ostream& in, Trio<T1, T2, T3>&);



}
}

#endif /* TRIO_HPP_ */
