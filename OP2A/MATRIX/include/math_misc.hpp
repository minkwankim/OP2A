/************************************************************************
			Open-source Multi-Physics Solver - ver. 0.0		
				
Copyright(c) 2013 MINKWAN KIM												
							
Initial Developed Date: Oct 29, 2013
                    by: Minkwan Kim

Last Modified Date: Oct 28, 2013
					by: Minkwan Kim


math_misc
	- Additional math functions 
              
					Copyright 2013 MINKWAN KIM
*************************************************************************/

#ifndef _MATH_MISC_H
#define _MATH_MISC_H

#include <math.h>
#include "../../UTIL/include/OP2A_utilities.hpp"

#include "matrix.hpp"
#include "vector_MK.hpp"



#ifndef MATH_PI
#define MATH_PI		3.14159265359
#endif


template <class T> T fabs(T a);

// AREA/LENGTH FUNCTIONS
double length(double n1[3], double n2[3], int ND);
double area_triangle(double n1[3], double n2[3], double n3[3], int ND);
double area_quadrilateral(double n1[3], double n2[3], double n3[3], double n4[3], int ND);
double volume_tetrahedron(double n1[3], double n2[3], double n3[3], double n4[3]);
double volume_pyramid(double n1[3], double n2[3], double n3[3], double n4[3], double n5[3]);
double volume_wedge(double n1[3], double n2[3], double n3[3], double n4[3], double n5[3], double n6[3]);
double volume_hexahedron(double n1[3], double n2[3], double n3[3], double n4[3], double n5[3], double n6[3], double n7[3], double n8[3]);

// Misc functions
double fsgn(double val);
double fmin(double a, double b);
double fmax(double a, double b);
double min_vec(double *data, int N);
double max_vec(double *data, int N);
double Delta_fn(int i, int j);
double interpolate (double *x, double *y, int n, double xp);
double interpolate2 (vector<double> &x, vector<double> &y, int n, double xp);
double MK_fn_ver1(int Xa, int Xb, double y_min, double y_max, int n);


#endif
