/************************************************************************
			Open-source Multi-Physics Solver - ver. 0.0		
				
Copyright(c) 2013 MINKWAN KIM												
							
Initial Developed Date: Oct 29, 2013
                    by: Minkwan Kim

Last Modified Date: Oct 28, 2013
					by: Minkwan Kim


CVector_MK.h
	- Library for CVector_MK calculations
              
					Copyright 2013 MINKWAN KIM
*************************************************************************/

#ifndef _CVector
#define _CVector

#include <math.h>

class CVector_MK
{

public:
	double xyz[3];	// Element of CVector_MK

	//////////////////////////////////
	/* Member function declarations */
	//////////////////////////////////
	// Constructor
	inline CVector_MK();
	inline CVector_MK(double x, double y, double z);
	inline CVector_MK(CVector_MK const &other);

	//::M-01
	//len()
	inline double len() const;

	//::M-02
	//normalize()
	inline CVector_MK normalize();

	//::M-03
	//normal_from _three_point
	inline CVector_MK normal_from_three_point(CVector_MK const &a, CVector_MK const &b, CVector_MK const &c);

	//::M-04
	// operator *=()
	inline void operator *= (double const & s);
	
	//::M-05
	// operator /=()
	inline void operator /= (double const & s);

	//::M-06a
	// operator -=()
	inline void operator += (CVector_MK const &s);

	//::M-06b
	// operator -=()
	inline void operator += (double const &s);

	//::M-07a
	// operator -=()
	inline void operator -= (CVector_MK const &s);
	
	//::M-07b
	// operator -=()
	inline void operator -= (double const &s);

	//::M-08
	// SetX()
	inline void SetX(double const & a);

	//::M-09
	// SetY()
	inline void SetY(double const & a);

	//::M-10
	// SetZ()
	inline void SetZ(double const & a);

	//::M-11
	// x()
	inline double x() const;

	//::M-12
	// y()
	inline double y() const;

	//::M-13
	// z()
	inline double z() const;

	//::M-14
	// Operator ()
	inline void operator()(double const x, double const y, double const z);
};




///////////////////////////////////
/* Global functions declarations */
///////////////////////////////////

//::01
//operator -()
inline CVector_MK operator -(CVector_MK const &a, CVector_MK const &b);

//::02
//operator %(): dot product
inline double operator %(CVector_MK const &a, CVector_MK const &b);

//::03
// operator *()
inline CVector_MK operator *(CVector_MK const &a, double const &s);

//::04
// operator *()
inline CVector_MK operator *(double const &s, CVector_MK const &a);

//::05
// operator *(): cross product
inline CVector_MK operator *(CVector_MK const &a, CVector_MK const &b);

//::06
// operator /() 
inline CVector_MK operator /(CVector_MK const & a, double const &s);

//::07
// operator +()
inline CVector_MK operator +(CVector_MK const & a, CVector_MK const &b);






//////////////////////////////////
/* Global function definitions  */
//////////////////////////////////

//::01
//operator -()
inline CVector_MK operator -(CVector_MK const & a, CVector_MK const &b)
{
	return(CVector_MK(a.x() - b.x(), a.y() - b.y(), a.z() - b.z()));
}

//::02
//operator %()
inline double operator %(CVector_MK const &a, CVector_MK const &b)
{
	return(a.x()*b.x() + a.y()*b.y() + a.z()*b.z());
}

//::03
// operator *()
inline CVector_MK operator *(CVector_MK const &a, double const &s)
{
	 return(CVector_MK(a.x()*s, a.y()*s, a.z()*s));
}

//::04
// operator *()
inline CVector_MK operator *(double const &s, CVector_MK const &a)
{
	return(CVector_MK(a.x()*s, a.y()*s, a.z()*s));
}

//::05
// operator *(): cross product
inline CVector_MK operator *(CVector_MK const &a, CVector_MK const &b)
{
	double  aa = a.y()*b.z()-a.z()*b.y();
	double  bb = b.x()*a.z()-a.x()*b.z();
	double  cc = a.x()*b.y()-a.y()*b.x();
	
	return(CVector_MK(aa, bb,cc));
}

//::06
// operator /() 
inline CVector_MK operator /(CVector_MK const & a, double const &s)
{
	return(CVector_MK(a.x()/s, a.y()/s, a.z()/s));
}

//::07
// operator +()
inline CVector_MK operator +(CVector_MK const & a, CVector_MK const &b)
{
	return(CVector_MK(a.x() + b.x(), a.y() + b.y(), a.z() + b.z() ));
}




/////////////////////////////////
/* Member function definitions */
/////////////////////////////////

// Constructor
inline CVector_MK::CVector_MK(double x, double y, double z)
{
	xyz[0] = x;
	xyz[1] = y;
	xyz[2] = z;
}

inline CVector_MK::CVector_MK(CVector_MK const &other)
{
	xyz[0] = other.x();
	xyz[1] = other.y();
	xyz[2] = other.z();
}

inline CVector_MK::CVector_MK()
{
	xyz[0] = 0.0;
	xyz[1] = 0.0;
	xyz[2] = 0.0;
}

//::M-01
//len()
inline double CVector_MK::len() const
{
	return( sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]) );
}

//::M-02
//normalize()
inline CVector_MK CVector_MK::normalize()
{
	 double len;
	 
	 len = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]) ;
	 xyz[0] /=len;
	 xyz[1] /=len;
	 xyz[2] /=len;
	 return  CVector_MK(xyz[0], xyz[1], xyz[2]);
}

//::M-03
//normal_from _three_point
inline CVector_MK CVector_MK::normal_from_three_point(CVector_MK const &a, CVector_MK const &b, CVector_MK const &c)
{
    CVector_MK ab=b-a;
	CVector_MK ac=c-a;
	return ab*ac;
}

//::M-04
// operator *=()
inline void CVector_MK::operator *=(double const & s)
{
	xyz[0] *= s;
	xyz[1] *= s;
	xyz[2] *= s;
}

//::M-05
// operator /=()
inline void CVector_MK::operator /=(double const & s)
{
	 xyz[0] /= s;
	 xyz[1] /= s;
	 xyz[2] /= s;
}

//::M-06a
// operator +=()
inline void CVector_MK::operator +=(CVector_MK const &s)
{
	xyz[0] += s.x();
	xyz[1] += s.y();
	xyz[2] += s.z();
}

//::M-06b
// operator +=()
inline void CVector_MK::operator +=(double const &s)
{
	 xyz[0] +=s;
	 xyz[1] +=s;
	 xyz[2] +=s;
}

//::M-07a
// operator -=()
inline void CVector_MK::operator -=(double const &s)
{
	 xyz[0] -=s;
	 xyz[1] -=s;
	 xyz[2] -=s;
}

//::M-07b
// operator -=()
inline void CVector_MK::operator -=(CVector_MK const &s)
{
	 xyz[0] -= s.x();
	 xyz[1] -= s.y();
	 xyz[2] -= s.z();
}

//::M-08
// SetX()
inline void CVector_MK::SetX(double const & a)
{
	 xyz[0]=a;
}

//::M-09
// SetY()
inline void CVector_MK::SetY(double const & a)
{
	 xyz[1]=a;
}

//::M-10
// SetZ()
inline void CVector_MK::SetZ(double const & a)
{
	 xyz[2]=a;
}

//::M-11
// x()
inline double CVector_MK::x() const
{
	return(xyz[0]);
}

//::M-12
// y()
inline double CVector_MK::y() const
{
	return(xyz[1]);
}

//::M-13
// z()
inline double CVector_MK::z() const
{
	return(xyz[2]);
}

//::M-14
// Operator ()
inline void CVector_MK::operator()(double const x, double const y, double const z)
{
	xyz[0] = x;
	xyz[1] = y;
	xyz[2] = z;
}

#endif 
