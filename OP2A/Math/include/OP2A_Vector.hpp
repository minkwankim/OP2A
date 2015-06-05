/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 5, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_Vector.hpp
 * 			-  
 *  
 */
#ifndef OP2A_VECTOR_HPP_
#define OP2A_VECTOR_HPP_


#include "Common/include/OP2A.hpp"

namespace OP2A{
namespace Math{


class VECTOR
{
public:
	VECTOR();
	VECTOR(const double x, const double y);
	VECTOR(const double x, const double y, double z);
	VECTOR(const std::vector<double>& x);
	VECTOR(const std::vector<double>& s, const std::vector<double>& e);

	~VECTOR();

public:	//Basic functions
	bool located();
	unsigned int size() const;
	unsigned int dimension() const;
	void resize(unsigned int s);
	double length();
	void normalize();
	double val(const unsigned int i) const;

public:	// Operators
	double& operator() (const unsigned int i);
	void operator *= (double const & s);
	void operator /= (double const & s);
	void operator -= (double const & s);
	void operator += (double const & s);

	void operator -= (const VECTOR & S);
	void operator += (const VECTOR & S);

public: // Linear Algebra



protected:
	unsigned int ND;
	std::vector<double>	m_data;
	bool is_allocated;


};

// Operators Functions
VECTOR operator- (const VECTOR &A, const VECTOR &B);
VECTOR operator+ (const VECTOR &A, const VECTOR &B);
VECTOR operator* (const double a, const VECTOR &B);



// Function for linear algebra
double VectorDotProduct(const VECTOR &A, const VECTOR &B);
VECTOR VectorCrossProduct(const VECTOR &A, const VECTOR &B);
VECTOR NormalFromThreePoint(const VECTOR &A, const VECTOR &B, const VECTOR &C);

}
}


#endif /* OP2A_VECTOR_HPP_ */
