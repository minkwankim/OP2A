/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 5, 2015
 *      			Author: Minkwan Kim
 *
 * MathMisc.hpp
 * 			-  
 *  
 */
#ifndef MATHMISC_HPP_
#define MATHMISC_HPP_

#include <vector>
#include "Common/include/OP2A.hpp"
//#include "Math/include/OP2A_Math.hpp"

namespace OP2A{
namespace Math{

#define MATH_PI		3.14159265359
#define MATH_ZERO	1.0e-15


// Misc functions
template <class T> T fabs(T a)
{
	if (a < 0) a = -a;
	return a;
}


template <typename T>
T fsgn(T val)
{
	T	plus	= 1;
	T	minus	= -1;

	if (val < 0)		return (minus);
	else if (val > 0)	return (plus);

    return (0);
}


template <typename T>
T fmin(T a, T b)
{
	if (a <= b)	return (a);
	else 		return (b);
}

template <typename T>
T fmin(std::vector<T> a)
{
	int N	= a.size();

	T min_temp = fmin<T>(a[0], a[1]);
	for (int i = 2; i <= N-1; i++)	min_temp = fmin<T>(min_temp, a[i]);

	return (min_temp);
}


template <typename T>
T fmax(T a, T b)
{
	if (a >= b)	return (a);
	else 		return (b);
}

template <typename T>
T fmax(std::vector<T> a)
{
	int N	= a.size();

	T max_temp = fmax<T>(a[0], a[1]);
	for (int i = 2; i <= N-1; i++)	max_temp = fmax<T>(max_temp, a[i]);

	return (max_temp);
}

template <typename T>
T Delta_fn(int i, int j)
{
	if (i == j) return (1);
	else		return (0);
}


double fminOMP(const std::vector<double> &V, const int iStart, const int iEnd);
double fmaxOMP(const std::vector<double> &V, const int iStart, const int iEnd);


double interpolate1(std::vector<double> &x, std::vector<double> &y, double xp);
double MK_fn_ver1(int Xa, int Xb, double y_min, double y_max, int n);

void Solve3x3eqn(std::vector< std::vector<double> >& A, std::vector<double> &b,  std::vector<double> &x);


}
}


#endif /* MATHMISC_HPP_ */
