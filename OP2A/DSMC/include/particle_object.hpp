/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: Dec 12, 2014
 *      			Author: Minkwan Kim
 *
 * particle_object.hpp
 * 			-  
 *  
 */

#ifndef PARTICLE_OBJECT_HPP_
#define PARTICLE_OBJECT_HPP_

#include <stddef.h>

#include "particle.hpp"

class particle_object
{
public:
	particle_type 	prop;		// Particle properties
	particle_object	*next;		// Pointer to next object

	/*
	 * Constructor and destructors
	 */
	particle_object();
	~particle_object();


	// Internal-functions
};


extern	particle_object	*unused_object;



// External-functions
// Ext_F-01: get_object(particle_object *a)
//			- Get an object from the list of unused object if it is available
//			- or allocate memory for a new object
// Ext_F-01: getobj(particle_object& b_unused)
//			- Get an object from the list of unused object if it is available
//			- or allocate memory for a new object
template <class T>
void GET(T *a, T *unused_list)
{
	if (unused_list != NULL)
	{
		a				= unused_list;
		unused_list		= unused_list->next;
		a->next			= NULL;
	}
	else
	{
		a		= new T;
		a->next	= NULL;
	}
}

// Ext_F-02: addobj(particle_object *a, particle_object *list)
//			- Add an object to the head of a list
template <class T>
void ADD(T *a, T *list)
{
	a->next = list;
	list	= a;
}




#endif /* PARTICLE_OBJECT_HPP_ */
