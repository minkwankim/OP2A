/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 14, 2015
 *      			Author: Minkwan Kim
 *
 * NamedObject.hpp
 * 			-  
 *  
 */
#ifndef NAMEDOBJECT_HPP_
#define NAMEDOBJECT_HPP_

#include "Common/include/Common.hpp"


namespace OP2A{
namespace Common{

class Common_API NamedObject
{
public:
	explicit	NamedObject(const std::string& name = std::string()) : m_name(name)
	{

	}

	virtual		~NamedObject()
	{

	}

	std::string	getName() const
	{
		return m_name;
	}

protected:
	void	setName(const std::string& name)
	{
		m_name	= name;
	}

private:
	std::string	m_name;	// Object name



};





}
}




#endif /* NAMEDOBJECT_HPP_ */
