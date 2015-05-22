/************************************************************************
			Open-source Multi-Physics Solver - ver. 0.0		
				
Copyright(c) 2013 MINKWAN KIM												
							
Initial Developed Date: Oct 29, 2013
                    by: Minkwan Kim

Last Modified Date: Oct 28, 2013
					by: Minkwan Kim


vector_2D.h
	- Library for vector calculations
              
					Copyright 2013 MINKWAN KIM
*************************************************************************/

#ifndef _VECTOR2D
#define _VECTOR2D

#include <math.h>

class vector2d
{
public:
	double xyz[2];	// Element of vector

	//////////////////////////////////
	/* Member function declarations */
	//////////////////////////////////
	// Constructor
	inline vector2d();
	inline vector2d(double a, double b);
	inline vector2d(vector2d const &other);

	//::M-01
	//len()
	inline double len() const;

	//::M-02
	//normalize()
	inline void normalize();

	//::M-03
	// x()/ y()
	inline double x() const;
	inline double y() const;

	//::M-04
	// SetX() / SetY()
	inline void SetX(double const & a);
	inline void SetY(double const & a);

	//::M-05
	// normal()
	inline vector2d normal(double a[2], double b[2]);
	inline void normal(vector2d const &A);

	//::M-06

};

///////////////////////////////////
/* Global functions declarations */
///////////////////////////////////
//::01
//operator -()
inline vector2d operator -(vector2d const &a, vector2d const &b);

//::02
//operator %(): dot product
inline double operator %(vector2d const &a, vector2d const &b);

//::03
// operator *()
inline vector2d operator *(vector2d const &a, double const &s);

//::04
// operator *()
inline vector2d operator *(double const &s, vector2d const &a);

//::05
// operator /() 
inline vector2d operator /(vector2d const & a, double const &s);

//::06
// operator +()
inline vector2d operator +(vector2d const & a, vector2d const &b);



//////////////////////////////////
/* Global function definitions  */
//////////////////////////////////

//::01
//operator -()
inline vector2d operator -(vector2d const & a, vector2d const &b)
{
	return(vector2d(a.x() - b.x(), a.y() - b.y()));
}

//::02
//operator %(): dot product
inline double operator %(vector2d const &a, vector2d const &b)
{
	return (a.x()*b.x() + a.y()*b.y());
}

//::03
// operator *()
inline vector2d operator *(vector2d const &a, double const &s)
{
	return (vector2d(a.x()*s, a.y()*s));
}

//::04
// operator *()
inline vector2d operator *(double const &s, vector2d const &a)
{
	return (vector2d(a.x()*s, a.y()*s));
}

//::05
// operator /() 
inline vector2d operator /(vector2d const & a, double const &s)
{
	return (vector2d(a.x()/s, a.y()/s));
}

//::06
// operator +()
inline vector2d operator +(vector2d const & a, vector2d const &b)
{
	return(vector2d(a.x()+b.x(), a.y()+b.y()));
}


/////////////////////////////////
/* Member function definitions */
/////////////////////////////////

// Constructor
inline vector2d::vector2d()
{
	xyz[0] = 0.0;
	xyz[1] = 0.0;
}
inline vector2d::vector2d(double a, double b)
{
	xyz[0] = a;
	xyz[1] = b;
}
inline vector2d::vector2d(vector2d const &other)
{
	xyz[0] = other.x();
	xyz[1] = other.y();
}



//::M-01
//len()
inline double vector2d::len() const
{
	return( sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]) );
}

//::M-02
//normalize()
inline void vector2d::normalize()
{
	 double length;
	 
	 length = len();
	 xyz[0] /=length;
	 xyz[1] /=length;
}

//::M-03
// x()
inline double vector2d::x() const
{
	return(xyz[0]);
}

// y()
inline double vector2d::y() const
{
	return(xyz[1]);
}

//::M-04
// SetX()
inline void vector2d::SetX(double const &a)
{
	xyz[0] = a;
}

// SetY()
inline void vector2d::SetY(double const &a)
{
	xyz[1] = a;
}

//::M-05
inline vector2d vector2d::normal(double a[2], double b[2])
{
	double length, Vx, Vy;
	Vx = b[0] - a[0];
	Vy = b[1] - a[1];

	length = sqrt(Vx*Vx + Vy*Vy);
	xyz[0] =  Vy / length;
	xyz[1] = -Vx / length;
	return vector2d(xyz[0], xyz[1]);
}

inline void vector2d::normal(vector2d const &A)
{
	double length;

	length = A.len();
	xyz[0] =  A.xyz[1] / length;
	xyz[1] = -A.xyz[0] / length;
}





#endif