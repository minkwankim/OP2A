/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 5, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_Vector.cpp
 * 			-  
 *  
 */



#include "Math/include/OP2A_Vector.hpp"
#include "Common/include/Exception_MemoryAllocator.hpp"
#include "Common/include/Exception_DimensionMatch.hpp"


namespace OP2A{
namespace Math{


VECTOR::VECTOR():ND(2), m_data(2, 0.0), is_allocated(false)
{

}

VECTOR::VECTOR(const unsigned int dim):ND(dim), m_data(dim, 0.0), is_allocated(true)
{

}


VECTOR::VECTOR(const double x, const double y):ND(2), m_data(2,0.0), is_allocated(true)
{
	m_data[0]	= x;
	m_data[1]	= y;
}


VECTOR::VECTOR(const double x, const double y,  const double z):ND(3), m_data(3,0.0), is_allocated(true)
{
	m_data[0]	= x;
	m_data[1]	= y;
	m_data[2]	= z;
}

VECTOR::VECTOR(const std::vector<double>& x):ND(x.size()), m_data(x), is_allocated(true)
{

}


VECTOR::VECTOR(const std::vector<double>& s, const std::vector<double>& e)
{
	if (s.size() != e.size())	throw Common::ExceptionDimensionMatch (FromHere(), "Dimenction of X and Y values do not match");

	ND	= e.size();
	m_data.resize(ND);

	for (int i = 0; i <= ND-1; i++)	m_data[i]	= e[i] - s[i];

	is_allocated = true;
}



VECTOR::~VECTOR()
{

}


////////////////////////////////////////////////////////////////////
bool VECTOR::located()
{
	return(is_allocated);
}


unsigned int VECTOR::size() const
{
	return ND;
}


unsigned int VECTOR::dimension() const
{
	return ND;
}


void VECTOR::resize(unsigned int s)
{
	m_data.resize(s, 0.0);
	ND	= s;
}


double VECTOR::length()
{
	double temp = 0.0;

	if (is_allocated == true)
	{
		for (int i = 0; i <= ND-1; i++)	temp += m_data[i]*m_data[i];
	}
	else
	{
		throw Common::ExceptionMemoryAllocator (FromHere(), "First, you need to allocate the size of VECTOR");
	}

	return (sqrt(temp));
}


void VECTOR::normalize()
{
	double len = length();

	for (int i = 0; i <= ND-1; i++)	m_data[i] /= len;
}


double VECTOR::val(const unsigned int i) const
{
	if (i > ND)	throw Common::ExceptionMemoryAllocator(FromHere(), "Exceed the data size of DataStorageVector");
	if (i == 0) throw Common::ExceptionMemoryAllocator(FromHere(), "Index should be non-zero integer");

	return m_data[i-1];
}





}
}
