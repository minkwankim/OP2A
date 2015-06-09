/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 28, 2015
 *      			Author: Minkwan Kim
 *
 * Vector2D.hpp
 * 			-  
 *  
 */
#ifndef VECTOR2D_HPP_
#define VECTOR2D_HPP_


#include <vector>
#include "Common/include/OP2A.hpp"


namespace OP2A{
namespace Common{


template <class TYPE>
class Vector2D
{
public:
	Vector2D() : m_M(0), m_N(0)
	{

	}

	Vector2D(unsigned int M, unsigned int N): m_M(M), m_N(N), m_vectorData(N*M)
	{

	}

	Vector2D(unsigned int M, unsigned int N, TYPE value): m_M(M), m_N(N), m_vectorData(N*M, value)
	{

	}


	Vector2D(const Vector2D& other ):
		m_M(other.m_M),
		m_N(other.m_N),
		m_vectorData(other.m_vectorData)
	{

	}

	~Vector2D()
	{

	}


	/*
	 * Overloading operators
	 */
	TYPE& operator()(const unsigned i, const unsigned j)
	{
		op_assert (i <= m_M);
		op_assert (j <= m_N);

		return m_vectorData[(i-1) + m_M*(j-1)];
	}

	Vector2D<TYPE>& operator=(const Vector2D<TYPE> A)
	{
		m_M	= A.m_M;
		m_N = A.m_N;
		m_vectorData = A.m_vectorData;

		return *this;
	}


	/*
	 * Basic functions
	 */
	void resize(unsigned m, unsigned n)
	{
		op_assert(m > 0);
		op_assert(n > 0);

		m_M = m;
		m_N = n;
		m_vectorData.resize(m_M*m_M);
	}

	void SetRowSize(unsigned m)
	{
		op_assert(m > 0);
		resize(m, m_N);
	}

	void SetColSize(unsigned n)
	{
		op_assert(n > 0);
		resize(n, m_N);
	}

	unsigned int size_row()
	{
		return m_M;
	}

	unsigned int size_col()
	{
		return m_N;
	}

	void reserve (unsigned int m, unsigned int n)
	{
		op_assert(m > 0);
		op_assert(n > 0);

		m_vectorData.resize(m*n);
	}

	void clear()
	{
		resize(0, 0);
	}


	void addCol(unsigned int a)
	{
		std::vector <TYPE>  vectorData_temp((m_M+1)*m_N, 0.0);

		int inew, jnew;
		for (int i = 1; i <= m_M; i++)
		{
			if (i < a)	inew	 = i;
			else		inew	 = i+1;

			for (int j = 1; j <= m_N; j++)
			{
				jnew = j;
				vectorData_temp[(inew-1) + (m_M+1)*(jnew-1)] = m_vectorData[IJ(i,j)];
			}
		}

		resize(m_M+1, m_N);
		m_vectorData = vectorData_temp;
	}

	void addCol()
	{
		addCol(m_M+1);
	}


	void addRow(unsigned int b)
	{
		std::vector <TYPE>  vectorData_temp(m_M*(m_N+1), 0.0);

		int inew, jnew;

		for (int j = 1; j <= m_N; j++)
		{
			if (j < b)	jnew	 = j;
			else		jnew	 = j + 1;

			for (int i = 1; i <= m_M; i++)
			{
				inew = i;
				vectorData_temp[(inew-1) + m_M*(jnew-1)] = m_vectorData[IJ(i,j)];
			}
		}

		resize(m_M, m_N+1);
		m_vectorData = vectorData_temp;
	}

	void addRow()
	{
		addRow(m_N);
	}


	void deleteCol(unsigned int a)
	{
		op_assert (a <= m_M);
		std::vector <TYPE>  vectorData_temp((m_M-1)*m_N, 0.0);

		int inew, jnew;
		for (int i = 1; i <= m_M; i++)
		{
			if (i != a)
			{
				if (i < a)	inew	 = i;
				if (i > a)	inew	 = i-1;

				for (int j = 1; j <= m_N; j++)
				{
					jnew = j;
					vectorData_temp[(inew-1) + (m_M-1)*(jnew-1)] = m_vectorData[IJ(i,j)];
				}
			}
		}

		resize(m_M-1, m_N);
		m_vectorData = vectorData_temp;
	}


	void deleteRow(unsigned int b)
	{
		std::vector <TYPE>  vectorData_temp(m_M*(m_N-1), 0.0);

		int inew, jnew;

		for (int j = 1; j <= m_N; j++)
		{
			if (j != b)
			{
				if (j < b)	jnew	 = j;
				if (j > b)	jnew	 = j - 1;

				for (int i = 1; i <= m_M; i++)
				{
					inew = i;
					vectorData_temp[(inew-1) + m_M*(jnew-1)] = m_vectorData[IJ(i,j)];
				}
			}
		}

		resize(m_M, m_N-1);
		m_vectorData = vectorData_temp;
	}







protected:
	unsigned int	m_M;
	unsigned int	m_N;
	std::vector <TYPE>  m_vectorData;

private:
	unsigned int IJ(unsigned int i, unsigned int j)
	{
		op_assert (i <= m_M);
		op_assert (j <= m_N);

		return ((i-1) + m_M*(j-1));
	}
};




}
}




#endif /* VECTOR2D_HPP_ */
