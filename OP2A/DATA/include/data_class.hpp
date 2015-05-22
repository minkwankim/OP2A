/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 12, 2015
 *      			Author: Minkwan Kim
 *
 * data_class.hpp
 * 			-  
 *  
 */

#ifndef DATA_CLASS_HPP_
#define DATA_CLASS_HPP_

#include <stddef.h>
#include <vector>
#include <string>

#include "../../UTIL/include/OP2A_utilities.hpp"

using namespace std;


/*
 * Basic Data class
 */
class DATA_BASIC_ver1
{
public:
	unsigned int	NV[4];			// Number of variables
	vector<double>	Q;				// Data series 1
	vector<double>	V;				// Data series 2
	vector<double>	W;				// Data series 3
	vector<double>	data;


	// Constructor /destructor
	DATA_BASIC_ver1();
	~DATA_BASIC_ver1();

	// Internal function
	void allocate_size();
};

class DATA_BASIC_ver2
{
public:
	int				ID;		// ID number
	unsigned int	NV;		// number of variable
	vector<double>	data;	// Data series 1

	// Constructor /destructor
	DATA_BASIC_ver2();
	~DATA_BASIC_ver2();

	// Internal function
	void allocate_size();
};


class DATA_BASIC_2D
{
public:
	int							ID;					// ID number
	vector < vector<double> >	data;	// Data series 1

	// Constructor /destructor
	DATA_BASIC_2D();
	~DATA_BASIC_2D();

	// Internal function
	void allocate_size(int n, int m);
};

#define DATA_BASIC	DATA_BASIC_ver2



class SOL_CLASS_DATA_2D
{
public:
	vector<int>				whereis;	// location of data 		(data pointer/Task ID)
	vector<DATA_BASIC_2D *>	data_ptr;

	SOL_CLASS_DATA_2D();
	~SOL_CLASS_DATA_2D();

	/*
	 * Internal functions
	 */
	void resize(unsigned int ncm, unsigned int n, int m);
};



/*
 * Solution Data class
 */
class SOL_CLASS_DATA
{
public:
	vector<int>				whereis;	// location of data 		(data pointer/Task ID)
	vector<DATA_BASIC *>	data_ptr;

	SOL_CLASS_DATA();
	~SOL_CLASS_DATA();

	/*
	 * Internal functions
	 */
	void resize(unsigned int ncm, unsigned int nv);
};




/*
 * Solution class
 */
class SOL_CLASS_BASIC: public SOL_CLASS_DATA
{
public:
	int	NDATA;							//	total number of data
	int NV;
	vector<string>	variable_names;

	SOL_CLASS_BASIC();
	~SOL_CLASS_BASIC();

	// International functions
	void allocate_size(unsigned int ncm, unsigned int nv);
};




#define SOL_CLASS	SOL_CLASS_BASIC


#endif /* DATA_CLASS_HPP_ */
