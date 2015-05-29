/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 11, 2015
 *      			Author: Minkwan Kim
 *
 * Exception.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_HPP_
#define EXCEPTION_HPP_


#include "Common/include/CodeLocation.hpp"
#include "Common/include/NonCopyable.hpp"



namespace OP2A{

namespace Common{

/* Exception Manager class
 *
 * @author	Minkwan Kim
 * @version 1.0 20/05/2015
 */
class Common_API ExceptionManager : public Common::NonCopyable<ExceptionManager>
{
public:
	ExceptionManager();

	static ExceptionManager& getInstance ();

	bool	ExceptionOutputs;
	bool	ExceptionDumps;
	bool	ExceptionAborts;
};



/* Exception class
 *
 * @author	Minkwan Kim
 * @version 1.0 20/05/2015
 */
class Common_API Exception : public std::exception
{
public:
	virtual ~Exception () throw ();

	void append (const std::string& add) throw ();		// Append additional description into m_what
	const std::string& str () const throw ();				// Get contents of description (m_what)
	const char* what () const throw (); 					// Get contents of description (m_what) [char]

	std::string full_description ()	const throw ();	// Show full description of exception

	virtual std::string getClassName () throw ()
	{
		return m_class_name;
	}

protected:
	Code_location	m_where;
	std::string		m_msg;
	std::string		m_class_name;
	std::string		m_what;

	Exception (Code_location where, std::string msg, std::string className) throw ();
};



}
}


#endif /* EXCEPTION_HPP_ */
