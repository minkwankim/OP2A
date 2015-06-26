/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 26, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_MatrixBasic.cpp
 * 			-  
 *  
 */

#include "Math/include/OP2A_Matrix.hpp"
#include "Common/include/Exception_MemoryAllocator.hpp"
#include "Common/include/Exception_DimensionMatch.hpp"

using namespace std;

namespace OP2A{
namespace Math{

/*
 * Constructors
 * @author Minkwan Kim
 * @version 1.0 26/6.2015
 */
MATRIX::MATRIX():m_allocated(false), m_matlabType(true)
{
	m_I = 0;
	m_J = 0;
}

MATRIX::MATRIX(bool matlabType):m_allocated(false), m_matlabType(matlabType)
{
	m_I = 0;
	m_J = 0;
}

MATRIX::MATRIX(const unsigned int I, const unsigned int J)
	:m_allocated(true), m_matlabType(true), m_I(I), m_J(J), m_data(I, vector<double>(J, 0.0))
{

}

MATRIX::MATRIX(const unsigned int I, const unsigned int J, bool matlabType)
	:m_allocated(true), m_matlabType(matlabType), m_I(I), m_J(J), m_data(I, vector<double>(J, 0.0))
{

}

MATRIX::MATRIX(const std::vector< std::vector< double> >& data_matrix)
	:m_allocated(true), m_matlabType(true), m_I(data_matrix.size()), m_J(data_matrix[0].size()), m_data(data_matrix)
{

}

MATRIX::MATRIX(const std::vector< std::vector< double> >& data_matrix, bool matlabType)
	:m_allocated(true), m_matlabType(matlabType), m_I(data_matrix.size()), m_J(data_matrix[0].size()), m_data(data_matrix)
{

}


MATRIX::~MATRIX()
{

}



/*
 * Basic functions
 * @author Minkwan Kim
 * @version 1.0 26/6.2015
 */
bool MATRIX::is_MatlabType() const
{
	return (m_matlabType);
}

bool MATRIX::is_allocated() const
{
	return (m_allocated);
}

unsigned int MATRIX::sizeI() const
{
	if (m_allocated != true)
	{
		throw Common::ExceptionMemoryAllocator (FromHere(), "First, you need to allocate the size of MATRIX");
	}

	return (m_I);
}

unsigned int MATRIX::sizeJ() const
{
	if (m_allocated != true)
	{
		throw Common::ExceptionMemoryAllocator (FromHere(), "First, you need to allocate the size of MATRIX");
	}

	return (m_J);
}


void MATRIX::resize(unsigned int I, unsigned int J)
{
	m_data.resize(I);
	for (int i = 0; i <= I-1; i++)	m_data[i].resize(J);

	m_I	= I;
	m_J	= J;
}

double MATRIX::val(unsigned int I, unsigned int J) const
{
	if (m_matlabType == true)
	{
		if (I > m_I)	throw Common::ExceptionMemoryAllocator(FromHere(), "Exceed the row size of a MATRIX");
		if (I == 0) 	throw Common::ExceptionMemoryAllocator(FromHere(), "Matlab-Type MATRIX Index should be non-zero integer");

		if (J > m_J)	throw Common::ExceptionMemoryAllocator(FromHere(), "Exceed the row size of a MATRIX");
		if (J == 0) 	throw Common::ExceptionMemoryAllocator(FromHere(), "Matlab-Type MATRIX Index should be non-zero integer");

		return (m_data[I-1][J-1]);
	}
	else
	{
		if (I >= m_I)	throw Common::ExceptionMemoryAllocator(FromHere(), "Exceed the row size of a C-type MATRIX");
		if (J >= m_J)	throw Common::ExceptionMemoryAllocator(FromHere(), "Exceed the row size of a C-type MATRIX");

		return (m_data[I][J]);
	}
}



void MATRIX::ones()
{
	if (m_allocated != true)
	{
		throw Common::ExceptionMemoryAllocator (FromHere(), "First, you need to allocate the size of MATRIX");
	}

#pragma omp parallel for
	for (int i = 0; i <= m_I-1; i++)
	{
#pragma ivdep
		for (int j = 0; j <= m_J-1; j++)
		{
			m_data[i][j] = 1.0;
		}
	}
}

void MATRIX::ones(unsigned int I, unsigned int J)
{
	resize(I, J);

#pragma omp parallel for
	for (int i = 0; i <= m_I-1; i++)
	{
#pragma ivdep
		for (int j = 0; j <= m_J-1; j++)
		{
			m_data[i][j] = 1.0;
		}
	}
}

void MATRIX::zeros()
{
	if (m_allocated != true)
	{
		throw Common::ExceptionMemoryAllocator (FromHere(), "First, you need to allocate the size of MATRIX");
	}

#pragma omp parallel for
	for (int i = 0; i <= m_I-1; i++)
	{
#pragma ivdep
		for (int j = 0; j <= m_J-1; j++)
		{
			m_data[i][j] = 0.0;
		}
	}
}

void MATRIX::zeros(unsigned int I, unsigned int J)
{
	resize(I, J);

#pragma omp parallel for
	for (int i = 0; i <= m_I-1; i++)
	{
#pragma ivdep
		for (int j = 0; j <= m_J-1; j++)
		{
			m_data[i][j] = 0.0;
		}
	}
}


void MATRIX::diag()
{
	if (m_allocated != true)
	{
		throw Common::ExceptionMemoryAllocator (FromHere(), "First, you need to allocate the size of MATRIX");
	}

	if (m_I != m_J)
	{
		throw Common::ExceptionMemoryAllocator (FromHere(), "You need a square matrix to create a diagonal matrix");
	}


#pragma omp parallel for
	for (int i = 0; i <= m_I-1; i++)
	{
		m_data[i][i] = 1.0;
	}
}


void MATRIX::diag(unsigned int I)
{
	resize(I, I);

#pragma omp parallel for
	for (int i = 0; i <= m_I-1; i++)
	{
		m_data[i][i] = 1.0;
	}
}

void MATRIX::diag(unsigned int I, double value)
{
	resize(I, I);

#pragma omp parallel for
	for (int i = 0; i <= m_I-1; i++)
	{
		m_data[i][i] = value;
	}
}



void MATRIX::removeRow(unsigned int I)
{
	if (m_allocated != true)	throw Common::ExceptionMemoryAllocator (FromHere(), "First, you need to allocate the size of MATRIX");

	if (m_matlabType == true)
	{
		if (I > m_I)	throw Common::ExceptionMemoryAllocator (FromHere(), "It is bigger than a current matrix size. Cannot remove a row");
		if (I == 0) 	throw Common::ExceptionMemoryAllocator(FromHere(), "removeRow: Matlab-Type MATRIX Index should be non-zero integer");

		int inew, jnew;
		vector <vector <double> >  data_temp(m_I-1, vector<double>(m_J)) ;

#pragma omp parallel for private(inew, jnew)
		for (int i = 0; i <= m_I-1; i++)
		{
			if (i != I-1)
			{
				if (i < I-1)	inew	 = i;
				if (i > I-1)	inew	 = i - 1;

#pragma ivdep
				for (int j = 0; j <= m_J-1; j++)
				{
					data_temp[inew][j] = m_data[i][j];
				}
			}
		}

		resize(m_I-1, m_J);
		m_data = data_temp;
	}
	else
	{
		if (I >= m_I)	throw Common::ExceptionMemoryAllocator (FromHere(), "It is bigger than a current matrix size. Cannot remove a row");

		int inew, jnew;
		vector <vector <double> >  data_temp(m_I-1, vector<double>(m_J)) ;

#pragma omp parallel for private(inew, jnew)
		for (int i = 0; i <= m_I-1; i++)
		{
			if (i != I)
			{
				if (i < I)	inew	 = i;
				if (i > I)	inew	 = i - 1;

#pragma ivdep
				for (int j = 0; j <= m_J-1; j++)
				{
					data_temp[inew][j] = m_data[i][j];
				}
			}
		}

		resize(m_I-1, m_J);
		m_data = data_temp;
	}
}



void MATRIX::removeColumn(unsigned int J)
{
	if (m_allocated != true)	throw Common::ExceptionMemoryAllocator (FromHere(), "First, you need to allocate the size of MATRIX");

	if (m_matlabType == true)
	{
		if (J > m_J)	throw Common::ExceptionMemoryAllocator (FromHere(), "It is bigger than a current matrix size. Cannot remove a column");
		if (J == 0) 	throw Common::ExceptionMemoryAllocator(FromHere(), "removeColumn: Matlab-Type MATRIX Index should be non-zero integer");

		int inew, jnew;
		vector <vector <double> >  data_temp(m_I, vector<double>(m_J-1)) ;

#pragma omp parallel for private(inew, jnew)
		for (int i = 0; i <= m_I-1; i++)
		{

#pragma ivdep
			for (int j = 0; j <= m_J-1; j++)
			{
				if (j != J-1)
				{
					if (j < J-1)	jnew	 = j;
					if (j > J-1)	jnew	 = j - 1;

					data_temp[i][jnew] = m_data[i][j];
				}
			}
		}

		resize(m_I, m_J-1);
		m_data = data_temp;
	}
	else
	{
		if (J >= m_J)	throw Common::ExceptionMemoryAllocator (FromHere(), "It is bigger than a current matrix size. Cannot remove a column");

		int inew, jnew;
		vector <vector <double> >  data_temp(m_I, vector<double>(m_J-1)) ;

#pragma omp parallel for private(inew, jnew)
		for (int i = 0; i <= m_I-1; i++)
		{
#pragma ivdep
			for (int j = 0; j <= m_J-1; j++)
			{
				if (j != J)
				{
					if (j < J)	jnew	 = j;
					if (j > J)	jnew	 = j - 1;

					data_temp[i][jnew] = m_data[i][j];
				}
			}
		}

		resize(m_I, m_J-1);
		m_data = data_temp;
	}

}


double MATRIX::sum()
{
	double o_sum = 0.0;
	vector<double> o_sumMP(m_I, 0.0);

#pragma omp parallel for
	for (int i = 0; i <= m_I-1; i++)
	{
		double sum_temp;
#pragma omp parallel for reduction(+: sum_temp)
		for (int j = 0; j <= m_J-1; j++)
		{
			sum_temp += m_data[i][j];
		}

		o_sumMP[i] = sum_temp;
	}


#pragma omp parallel for reduction(+: o_sum)
	for (int i = 0; i <= m_I-1; i++)
	{
		o_sum += o_sumMP[i];
	}

	return (o_sum);
}




}
}
