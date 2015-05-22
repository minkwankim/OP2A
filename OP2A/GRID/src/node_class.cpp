/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 7, 2015
 *      			Author: Minkwan Kim
 *
 * node_class.cpp
 * 			-  
 *  
 */


#include "../include/node_class.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"

/*
 * Basic information of node class
 */
NODE_BASIC_2D::NODE_BASIC_2D()
{
	ID		= 0;
	x[0]	= 0.0;
	x[1]	= 0.0;
}

NODE_BASIC_2D::~NODE_BASIC_2D()
{

}

NODE_BASIC_3D::NODE_BASIC_3D()
{
	ID		= 0;
	x[0]	= 0.0;
	x[1]	= 0.0;
	x[2]	= 0.0;
}

NODE_BASIC_3D::~NODE_BASIC_3D()
{

}


/*
 * Extended information of node class
 */
NODE_OPTIONAL::NODE_OPTIONAL()
{
	NSC	= 0;
	C_list.reserve(4);
	W_C.reserve(4);
}

NODE_OPTIONAL::~NODE_OPTIONAL()
{

}

// F-01:: Assign asign of C_list
void NODE_OPTIONAL::assign_C_list_size()
{
	if (NSC > 0)
	{
		C_list.resize(NSC);
		W_C.resize(NSC);
	}
	else
	{
		Error_message_type	error;
		error.message		= "Number of shared cell(NSC) is not assigned yet.";
		error.module_name	= "node_class.cpp: NODE_OPTIONAL";
		error.print_message();
	}
}





/*
 * Node Class
 */
NODE_CLASS::NODE_CLASS()
{

}

NODE_CLASS::~NODE_CLASS()
{

}
