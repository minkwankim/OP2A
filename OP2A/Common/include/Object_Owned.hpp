/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * Object_Owned.hpp
 * 			-  
 *  
 */
#ifndef OBJECT_OWNED_HPP_
#define OBJECT_OWNED_HPP_

#include "Common/include/Common.hpp"
#include "Common/include/OP2A.hpp"


namespace OP2A{
namespace Common{

class Common_API ObjectOwned{
public:

	ObjectOwned() : m_owners(0) {}
	virtual ~ObjectOwned() {}

	/// Add an owner.
	void addOwner()	{ m_owners++; }

	/// Remove an owner.
	void removeOwner() { m_owners--; }

	/// Check if the object is not owned anymore.
	bool hasNoOwner() { return (m_owners == 0) ? true : false; }

	/// Get the number of owners
	OPuint getNbOwners() const { return m_owners; }

private:
	OPuint m_owners;
};

}// end namespace Common
}// end namespace OP2A



#endif /* OBJECT_OWNED_HPP_ */
