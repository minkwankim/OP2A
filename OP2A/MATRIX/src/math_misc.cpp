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
#include "../include/math_misc.hpp"


template <class T>
T fabs(T a)
{
	if (a < 0) a = -a;

	return a;
}

double fsgn(double val)
{
	double aux;
	if (val < 0.0)
	{
		aux	= -1.0;
	}
	else if (val > 0.0)
	{
		aux	= 1.0;
	}
	else
	{
		aux	= 0.0;
	}

    return (aux);
}

double fmin(double a, double b)
{
	if (a <= b)
	{
		return (a);
	}
	else
	{
		return (b);
	}
}

double fmax(double a, double b)
{
	if (a >= b)
	{
		return (a);
	}
	else
	{
		return (b);
	}
}


double min_vec(double *data, int N)
{
	double min	= data[0];
	for (int i = 1; i <= N-1; i++)
	{
		if (data[i] < min)	min = data[i];
	}

	return (min);
}

double max_vec(double *data, int N)
{
	double max	= data[0];
	for (int i = 1; i <= N-1; i++)
	{
		if (data[i] > max)	max = data[i];
	}

	return (max);
}


////////////////////////////////
/* AREA CALCULATION FUNCTIONS */
////////////////////////////////


double length(double n1[3], double n2[3], int ND)
{
	int i;
	double len;

	len = 0.0;
	for (i = 0; i <= ND-1; i++)
	{
		len += pow((n1[i] - n2[i]), 2.0);
	}
	return (sqrt(len));
}

double area_triangle(double n1[3], double n2[3], double n3[3], int ND)
{
	double a, b, c, s;
	double area;
    
	a = length(n1, n2, ND);
	b = length(n2, n3, ND);
	c = length(n3, n1, ND);

	s = 0.5 * (a + b + c);
	area = sqrt(s * (s-a) * (s-b) * (s-c));

	return (area);
}

double area_quadrilateral(double n1[3], double n2[3], double n3[3], double n4[3], int ND)
{
	double s1, s2;
	
	s1 = area_triangle(n1, n2, n3, ND);
	s2 = area_triangle(n2, n3, n4, ND);

	return (s1 + s2);
}


double volume_tetrahedron(double n1[3], double n2[3], double n3[3], double n4[3])
{
	int error_code;
	string error_message;

	double vol;
	
	CVector_MK a, b, c, d;
	CVector_MK A, B, C, TEMP1;

	a.SetX(n1[0]);
	a.SetY(n1[1]);
	a.SetZ(n1[2]);

	b.SetX(n1[0]);
	b.SetY(n1[1]);
	b.SetZ(n1[2]);

	c.SetX(n1[0]);
	c.SetY(n1[1]);
	c.SetZ(n1[2]);

	d.SetX(n1[0]);
	d.SetY(n1[1]);
	d.SetZ(n1[2]);

	A = a - d;
	B = b - d;
	C = c - d;

	TEMP1	= B * C;
	vol		= A % TEMP1;
	vol		= fabs(vol) / 6.0;

	if (error_check_double_pos(vol))
	{
		Error_message_type	error;
		error.location_primary_name	="math_misc.cpp";
		error.message = "PROBLEM IN VOLUME CALCULATION OF TETRAHEDRON!";
		error.print_message();
	}
	
	return (vol);
}

/* [NOTE] n5 is the top point of a pyramid */
double volume_pyramid(double n1[3], double n2[3], double n3[3], double n4[3], double n5[3])
{
	int i;
	double volume;
	double n6[3];

	int error_code;
	string error_message;

	for (i = 0; i <= 2; i++)
	{
		n6[i] = 0.25 * (n1[i] + n2[i] + n3[i] + n4[i]);
	}

	volume = volume_tetrahedron(n1, n2, n6, n5);
	volume += volume_tetrahedron(n2, n3, n6, n5);
	volume += volume_tetrahedron(n3, n4, n6, n5);
	volume += volume_tetrahedron(n4, n1, n6, n5);


	if (error_check_double_pos(volume))
	{
		Error_message_type	error;
		error.location_primary_name	="math_misc.cpp";
		error.message = "PROBLEM IN VOLUME CALCULATION OF PYRAMID!";
		error.print_message();
	}
	return (volume);
}

double volume_wedge(double n1[3], double n2[3], double n3[3], double n4[3], double n5[3], double n6[3])
{
	int i;
	double n7[3];
	double volume;

	int error_code;
	string error_message;

	for (i = 0; i <= 2; i++)
	{
		n7[i] = (n1[i] + n2[i] + n3[i] + n4[i] + n5[i] + n6[i]) / 6.0;
	}
	
	volume = volume_pyramid(n1, n4, n5, n2, n7);
	volume = volume + volume_pyramid(n2, n3, n6, n3, n7);
	volume = volume + volume_pyramid(n1, n3, n6, n4, n7);
	volume = volume + volume_tetrahedron(n1, n2, n3, n7);
	volume = volume + volume_tetrahedron(n4, n6, n5, n7);

	if (error_check_double_pos(volume))
	{
		Error_message_type	error;
		error.location_primary_name	="math_misc.cpp";
		error.message = "PROBLEM IN VOLUME CALCULATION OF WEDGE!";
		error.print_message();
	}
	return (volume);
}

double volume_hexahedron(double n1[3], double n2[3], double n3[3], double n4[3], double n5[3], double n6[3], double n7[3], double n8[3])
{
	int i;
	double n9[3];
	double volume;

	int error_code;
	string error_message;

	for (i = 0; i <= 2; i++)
	{
		n9[i] = 0.125 * (n1[i] + n2[i] + n3[i] + n4[i] + n5[i] + n6[i] + n7[i] + n8[i]);
	}

	volume = volume_pyramid(n1, n2, n3, n4, n9);
	volume += volume_pyramid(n5, n8, n7, n6, n9);
	volume += volume_pyramid(n1, n5, n6, n2, n9);
	volume += volume_pyramid(n3, n7, n8, n4, n9);
	volume += volume_pyramid(n2, n6, n7, n3, n9);
	volume += volume_pyramid(n1, n4, n8, n5, n9);

	if (error_check_double_pos(volume))
	{
		Error_message_type	error;
		error.location_primary_name	="math_misc.cpp";
		error.message = "PROBLEM IN VOLUME CALCULATION OF HEXAHEDRON!";
		error.print_message();
	}
	return (volume);
}





double Delta_fn(int i, int j)
{
	double a;
	if (i == j)	a = 1.0;
	else	a = 0.0;

	return a;
}


double interpolate (double *x, double *y, int n, double xp)
{
	int i;
	double x1,y1,x2,y2,yp;

	if (xp <= x[0])
	{
		yp=y[0];
	}
	else if (xp >= x[n-1])
	{
		yp=y[n-1];
	}
	else
	{
		for (i = 0; i <= n-2; i++)
		{
			if (xp >= x[i] && xp < x[i+1])
			{
				x1=x[i];
				y1=y[i];
				x2=x[i+1];
				y2=y[i+1];

				yp=y1+(y2-y1)/(x2-x1)*(xp-x1);

			}


		}


	}

	return(yp);
}


double interpolate2 (vector<double> &x, vector<double> &y, int n, double xp)
{
	int i;
	double x1,y1,x2,y2,yp;

	if (xp <= x[0])
	{
		yp=y[0];
	}
	else if (xp >= x[n-1])
	{
		yp=y[n-1];
	}
	else
	{
		for (i = 0; i <= n-2; i++)
		{
			if (xp >= x[i] && xp < x[i+1])
			{
				x1=x[i];
				y1=y[i];
				x2=x[i+1];
				y2=y[i+1];

				yp=y1+(y2-y1)/(x2-x1)*(xp-x1);

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





/*



inline int p_solve3x3eqn(double const a0, double const b0, double const c0, double const s0, 
				  double const a1, double const b1, double const c1, double const s1,
				  double const a2, double const b2, double const c2, double const s2,
				  double &u,   double &v,  double &w)
{	
    if( fabs(s0) < 1e-15   && fabs(s1) <1e-15 && fabs(s2) <1e-15) 
	{
	         u=v=w=0.0; 
			 return  1;
	}
	double det = a0*b1*c2 + a1*b2*c0 +a2*c1*b0 -a2*b1*c0 -a1*b0*c2 -a0*b2*c1; 	 
	if( fabs(det) <1e-17) 
		return -1; 
	
	double a11 = b1*c2-b2*c1; 
	double a12 = -(b0*c2-b2*c0); 
	double a13 = b0*c1-b1*c0; 

	double a21 =-(a1*c2-a2*c1); 
	double a22 = a0*c2-a2*c0; 
	double a23 =-(a0*c1-a1*c0); 

	double a31 = a1*b2-a2*b1; 
	double a32 =-(a0*b2-a2*b0); 
	double a33 = a0*b1-a1*b0; 

	u =  (a11*s0+a12*s1+a13*s2)/det; 
	v =  (a21*s0+a22*s1+a23*s2)/det;
	w =  (a31*s0+a32*s1+a33*s2)/det;

	return 1; 
}

inline int p_CalculateGradientbyLeastSquare(int const sizeN, double *dx, double *dy, double *dz, double *du)
{
	if( sizeN ==0) 
	{
		du[0]=du[1]=du[2]=0.0;
		return 1;
	}
	if( sizeN <3) 
		return -1; 

	double a0 =0., b0=0., c0=0., a1=0., b1=0., c1=0., a2=0., b2=0., c2=0.; 
	double rhs1=0., rhs2=0., rhs3=0.,  su0=0.0, su1=0.0, su2=0.0; 
	int  ii=0;
    for( ii=0; ii<sizeN; ii++) 
	{
		a0 += dx[ii] *dx[ii]; 
		b0 += dx[ii] *dy[ii]; 
		c0 += dx[ii] *dz[ii]; 

		a1 += dy[ii] *dx[ii]; 
		b1 += dy[ii] *dy[ii]; 
		c1 += dy[ii] *dz[ii]; 
		
		a2 += dz[ii] *dx[ii]; 
		b2 += dz[ii] *dy[ii]; 
		c2 += dz[ii] *dz[ii]; 

		rhs1 +=dx[ii] * du[ii]; 
        rhs2 +=dy[ii] * du[ii]; 
        rhs3 +=dz[ii] * du[ii]; 
	}

	double lmax = fabs(a0); 
	lmax = lmax>fabs(b0)? lmax:fabs(b0); 
	lmax = lmax>fabs(c0)? lmax:fabs(c0); 

	lmax = lmax>fabs(a1)? lmax:fabs(a1); 
	lmax = lmax>fabs(b1)? lmax:fabs(b1); 
	lmax = lmax>fabs(c1)? lmax:fabs(c1); 

	lmax = lmax>fabs(a2)? lmax:fabs(a2); 
	lmax = lmax>fabs(b2)? lmax:fabs(b2); 
	lmax = lmax>fabs(c2)? lmax:fabs(c2); 

	if( lmax <1e-3) 
	{
		lmax =1./lmax; 
		a0 *=lmax; 
		b0 *=lmax;
		c0 *=lmax; 

		a1 *=lmax; 
		b1 *=lmax; 
		c1 *=lmax; 

		a2 *=lmax; 
		b2 *=lmax; 
		c2 *=lmax; 

	    rhs1 *=lmax; 
		rhs2 *=lmax; 
		rhs3 *=lmax; 
	}

	int  ress =1;
	if( fabs(rhs1 <1e-20) && fabs(rhs2)<1e-20 && fabs(rhs3)<1e-20) 
		du[0]=du[1]=du[2]=0.0; 
	else 
		ress=p_solve3x3eqn(a0, b0, c0, rhs1, a1, b1, c1, rhs2, a2, b2, c2, rhs3, su0, su1, su2);

	if( ress <0) 
		return  ress; 

	du[0] =su0; 
	du[1] =su1; 
	du[2] =su2; 

	return 1; 
}

inline int  fourfourmatch(int n0, int n1, int n2, int n3,  int m0, int m1, int m2, int m3) 
{
	if(n0 == m0 || n0 ==m1 || n0 ==m2 || n0 ==m3) 
	{
			if(n1 == m0 || n1 ==m1 || n1 ==m2 || n1 ==m3)
			{
					if(n2 == m0 || n2 ==m1 || n2 ==m2 || n2 ==m3)
					{
						if(n3 == m0 || n3 ==m1 || n3 ==m2 || n3 ==m3)	 return 1; 
					}
			}
	}
	return 0; 
}
*/
