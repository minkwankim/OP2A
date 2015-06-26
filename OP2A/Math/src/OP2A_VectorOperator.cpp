/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 5, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_VectorOperator.hpp
 * 			-  
 *  
 */




#include "Math/include/OP2A_Vector.hpp"
#include "Common/include/Exception_MemoryAllocator.hpp"
#include "Common/include/Exception_DimensionMatch.hpp"


namespace OP2A{
namespace Math{

double& VECTOR::operator() (const unsigned int i)
{
	if (i > ND)	throw Common::ExceptionMemoryAllocator(FromHere(), "Exceed the data size of DataStorageVector");
	if (i == 0) throw Common::ExceptionMemoryAllocator(FromHere(), "Index should be non-zero integer");

	return m_data[i-1];
}


void VECTOR::operator *= (double const & s)
{
	for (int i = 0; i <= ND-1; i++)	m_data[i]	*= s;
}

void VECTOR::operator /= (double const & s)
{
	for (int i = 0; i <= ND-1; i++)	m_data[i]	/= s;
}

void VECTOR::operator += (double const & s)
{
	for (int i = 0; i <= ND-1; i++)	m_data[i]	+= s;
}

void VECTOR::operator -= (double const & s)
{
	for (int i = 0; i <= ND-1; i++)	m_data[i]	-= s;
}


void VECTOR::operator -= (const VECTOR& S)
{
	if (ND != S.ND)	throw Common::ExceptionDimensionMatch (FromHere(), "Dimension of X and Y values do not match");

	for (int i = 0; i <= ND-1; i++)	m_data[i]	-= S.m_data[i];
}

void VECTOR::operator += (const VECTOR& S)
{
	if (ND != S.ND)	throw Common::ExceptionDimensionMatch (FromHere(), "Dimension of X and Y values do not match");

	for (int i = 0; i <= ND-1; i++)	m_data[i]	+= S.m_data[i];
}





// Operator Functions
VECTOR operator- (const VECTOR &A, const VECTOR &B)
{
	unsigned int N = A.size();
	if (N != B.size())	throw Common::ExceptionDimensionMatch (FromHere(), "Dimension of X and Y values do not match");
	VECTOR C(N, 0.0);


	for (int i = 1; i <= N; i++)	C(i) = A.val(i) - B.val(i);
	return (C);
}


VECTOR operator+ (const VECTOR &A, const VECTOR &B)
{
	unsigned int N = A.size();
	if (N != B.size())	throw Common::ExceptionDimensionMatch (FromHere(), "Dimension of X and Y values do not match");
	VECTOR C(N, 0.0);


	for (int i = 1; i <= N; i++)	C(i) = A.val(i) + B.val(i);
	return (C);
}

VECTOR operator* (const double a, const VECTOR &B)
{
	unsigned int N = B.size();

	VECTOR C(N, 0.0);
	for (int i = 1; i <= N; i++)	C(i) = a * B.val(i);
	return (C);
}


}
}
