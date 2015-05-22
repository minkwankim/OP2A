/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 7, 2015
 *      			Author: Minkwan Kim
 *
 * node_class.hpp
 * 			-  
 *  
 */

#ifndef NODE_CLASS_HPP_
#define NODE_CLASS_HPP_

#include <vector>
using namespace std;

/*
 * Basic information of node class
 */
class NODE_BASIC_2D
{
public:
	int		ID;		// Node ID
	double x[2];

	// Constructor and destructor
	NODE_BASIC_2D();
	~NODE_BASIC_2D();
};

class NODE_BASIC_3D
{
public:
	int		ID;		// Node ID
	double x[3];

	// Constructor and destructor
	NODE_BASIC_3D();
	~NODE_BASIC_3D();
};

/*
 * Extended information of node class
 */
class NODE_OPTIONAL
{
public:
	unsigned int 	NSC;		// number of shared cells
	vector<double>	C_list;		// List of shared cells
	vector<double>	W_C;		// Weight of shared cells

	// Constructor and destructor
	 NODE_OPTIONAL();
	~ NODE_OPTIONAL();

	// F-01:: Assign asign of C_list
	void assign_C_list_size();
};


#define NODE_BASIC		NODE_BASIC_2D
#ifdef DIM_3D
#undef	NODE_BASIC
#define NODE_BASIC	NODE_BASIC_3D
#endif

/*
 * Node Class
 */
class NODE_CLASS : public NODE_BASIC, public NODE_OPTIONAL
{
public:
	NODE_CLASS();
	~NODE_CLASS();
};



#endif /* NODE_CLASS_HPP_ */
