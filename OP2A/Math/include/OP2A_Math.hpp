/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 5, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_Math.hpp
 * 			-  
 *  
 */
#ifndef OP2A_MATH_HPP_
#define OP2A_MATH_HPP_


#include "Common/include/OP2A.hpp"

namespace OP2A{
namespace Math{


#define MATH_PI		3.14159265359
#define MATH_ZERO	1.0e-15




template <class T> T fabs(T a)
{
	if (a < 0) a = -a;
	return a;
}




// AREA/LENGTH FUNCTIONS
double length(double n1[3], double n2[3], int ND);
double area_triangle(double n1[3], double n2[3], double n3[3], int ND);
double area_quadrilateral(double n1[3], double n2[3], double n3[3], double n4[3], int ND);
double volume_tetrahedron(double n1[3], double n2[3], double n3[3], double n4[3]);
double volume_pyramid(double n1[3], double n2[3], double n3[3], double n4[3], double n5[3]);
double volume_wedge(double n1[3], double n2[3], double n3[3], double n4[3], double n5[3], double n6[3]);
double volume_hexahedron(double n1[3], double n2[3], double n3[3], double n4[3], double n5[3], double n6[3], double n7[3], double n8[3]);


}
}

#endif /* OP2A_MATH_HPP_ */
