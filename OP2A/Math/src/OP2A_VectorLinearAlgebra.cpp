/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 5, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_VectorLinearAlgebra.cpp
 * 			-  
 *  
 */



#include "Math/include/OP2A_Vector.hpp"
#include "Common/include/Exception_MemoryAllocator.hpp"
#include "Common/include/Exception_DimensionMatch.hpp"


namespace OP2A{
namespace Math{


VECTOR VectorCrossProduct(const VECTOR &A, const VECTOR &B)
{
	unsigned int N = A.size();
	if (N != B.size() || N != 3)	throw Common::ExceptionDimensionMatch (FromHere(), "Dimenction of X and Y values do not match");
	VECTOR C(N, 0.0);

	C(1)	= A.val(2)*B.val(3)- A.val(3)*B.val(2);
	C(2)	= A.val(3)*B.val(1)- A.val(1)*B.val(3);
	C(3)	= A.val(1)*B.val(2)- A.val(2)*B.val(1);

	return (C);
}

double VectorDotProduct(const VECTOR &A, const VECTOR &B)
{
	unsigned int N = A.size();
	if (N != B.size())	throw Common::ExceptionDimensionMatch (FromHere(), "Dimenction of X and Y values do not match");

	double ans = 0.0;
	for (int i = 1; i <= N; i++)	ans	+= A.val(i)*B.val(i);
	return (ans);
}

VECTOR NormalFromThreePoint(const VECTOR &A, const VECTOR &B, const VECTOR &C)
{
	VECTOR AB = B - A;
	VECTOR AC = C - A;

	VECTOR ans;
	ans	= VectorCrossProduct(AB, AC);
	return (ans);
}


}
}
