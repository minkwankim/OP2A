/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 5, 2015
 *      			Author: Minkwan Kim
 *
 * MathMiscFn.hpp
 * 			-  
 *  
 */

#include <limits>
#include <algorithm>
#include "Common/include/Exception_DimensionMatch.hpp"
#include "Math/include/MathMisc.hpp"


namespace OP2A{
namespace Math{

double fminOMP(const std::vector<double> &V, const int iStart, const int iEnd)
{

	double Vmin = V[iStart];

#pragma omp parallel for reduction (min:Vmin)
	for (int i = iStart+1; i <= iEnd; i++)
	{
		if (V[i] < Vmin)	Vmin = V[i];
	}

	return Vmin;
}

double fmaxOMP(const std::vector<double> &V, const int iStart, const int iEnd)
{

	double Vmax = V[iStart];

#pragma omp parallel for reduction (min:Vmax)
	for (int i = iStart+1; i <= iEnd; i++)
	{
		if (V[i] > Vmax)	Vmax = V[i];
	}

	return Vmax;
}




double interpolate1(std::vector<double> &x, std::vector<double> &y, double xp)
{

	int i;
	int N	= x.size();

	if (N != y.size())	throw Common::ExceptionDimensionMatch (FromHere(), "Dimenction of X and Y values do not match");


	double x1, y1;
	double x2, y2;
	double yp;


	if (xp <= x[0])
	{
		yp = y[0];
	}
	else if (xp >= x[N-1])
	{
		yp = y[N-1];
	}
	else
	{
		for (i = 0; i <= N-2; i++)
		{
			if (xp >= x[i] && xp < x[i+1])
			{
				x1 = x[i];
				y1 = y[i];
				x2 = x[i+1];
				y2 = y[i+1];

				yp = y1 + (y2-y1)/(x2-x1)*(xp-x1);
			}
		}
	}

	return(yp);
}


double MK_fn_ver1(int Xa, int Xb, double y_min, double y_max, int n)
{
	double y;

	double temp1	= pow(Xa, 3.0)/3.0	- (Xa + Xb)/2.0*pow(Xa, 2.0)	+ pow(Xa, 2.0)*Xb;
	double temp2	= pow(Xb, 3.0)/3.0	- (Xa + Xb)/2.0*pow(Xb, 2.0)	+ pow(Xb, 2.0)*Xa;

	double A	= (y_min - y_max) / (-temp1 + temp2);
	double C	= y_max	+ temp2*A;
	double temp3	= Xa * Xb;

	if (n <= Xa)
	{
		y	= y_min;
	}
	else if (n >= Xb)
	{
		y	= y_max;
	}
	else
	{
		y	= -A* (pow(n, 3.0)/3.0	- (Xa + Xb)/2.0*pow(n, 2.0)	+ temp3*n) + C;
	}

	return (y);
}

void Solve3x3eqn(std::vector< std::vector<double> >& A, std::vector<double> &b,  std::vector<double> &x)
{
	if (fabs<double>(b[0]) <= MATH_ZERO	&& fabs<double>(b[1]) <= MATH_ZERO && fabs<double>(b[2]) <= MATH_ZERO)
	{
		x.resize(3, 0.0);
	}
	else
	{
		x.resize(3, 0.0);
		double det = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]) - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][1]) + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);

		if( fabs<double>(det) < MATH_ZERO) throw Common::ExceptionDimensionMatch (FromHere(), "NO Solution (zero determinent)");

		double a11 = A[1][1]*A[2][2] - A[1][2]*A[2][1];
		double a12 = -(A[1][0]*A[2][2] - A[1][2]*A[2][0]);
		double a13 = A[1][0]*A[2][1] - A[1][1]*A[2][0];

		double a21 =-(A[0][1]*A[2][2] - A[0][2]*A[2][1]);
		double a22 = A[0][0]*A[2][2] - A[0][2]*A[2][0];
		double a23 =-(A[0][0]*A[2][1] - A[0][1]*A[2][0]);

		double a31 = A[2][1]*A[1][2] - A[0][2]*A[1][1];
		double a32 =-(A[0][0]*A[1][2] - A[0][2]*A[1][0]);
		double a33 = A[0][0]*A[1][1] - A[0][1]*A[1][0];

		x[0] =  (a11*b[0] + a12*b[1] + a13*b[2])/det;
		x[1] =  (a21*b[0] + a22*b[1] + a23*b[2])/det;
		x[2] =  (a31*b[0] + a32*b[1] + a33*b[2])/det;
	}
}




}
}
