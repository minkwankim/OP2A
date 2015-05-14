/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 11, 2015
 *      			Author: Minkwan Kim
 *
 * Exception.cpp
 * 			-  
 *  
 */



#include "Common/include/Exception.hpp"
#include "Common/include/OSystem.hpp"
#include "Common/include/ProcessInfo.hpp"




namespace OP2A{

namespace Common{

ExceptionManager::ExceptionManager()
{
	ExceptionOutputs	= true;
	ExceptionDumps		= true;
	ExceptionAborts		= false;
}

ExceptionManager& ExceptionManager::getInstance()
{
	static	ExceptionManager exception_manager;
	return 	exception_manager;
}





Exception::Exception (Code_location where, std::string msg, std::string className) throw ()
		: m_where(where), m_msg(msg), m_class_name(className)
{
	m_what = full_description();

	if ( ExceptionManager::getInstance().ExceptionOutputs )
	{
		std::cout << "+++ Exception thrown +++++++++++++++++" << std::endl;
		std::cout << m_what << std::endl;
		std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
	}

	if ( ExceptionManager::getInstance().ExceptionDumps )
	{
		std::string backtrace = OSystem::getInstance().getProcessInfo()->getBackTrace();
		std::cout << "+++ Exception backtrace ++++++++++++++" << std::endl;
		std::cout << backtrace << std::endl;
		std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
	}

	if ( ExceptionManager::getInstance().ExceptionAborts )
	{
		std::cout << "+++ Exception aborting ... " << std::endl;
		abort();
	}
}


Exception::~Exception() throw ()
{

}

void Exception::append(const std::string& add) throw ()
{
	m_what += add;
}

const std::string& Exception::str () const throw ()
{
	return m_what;
}

const char* Exception::what() const throw()
{
	return str().c_str();
}

std::string Exception::full_description () const throw ()
{
	std::string desc;
	desc += "+++ Exception thrown ++++++++++++++++\n";
	desc += "From : \'";
	desc += m_where.str();
	desc += "\'\n";
	desc += "Type : \'";
	desc += m_class_name;
	desc += "\'\n";
	desc += "Message :\n";
	desc += m_msg;
	desc += "\n+++++++++++++++++++++++++++++++++++++\n";
	return desc;
}


std::ostream& operator<<(std::ostream& output, OP2A::Common::Exception& exc)
{
	output << exc.str();
	return output;
}

std::ostream& operator<<(std::ostream& output, OP2A::Common::Exception* exc)
{
  return output << *exc;
}


}
}
