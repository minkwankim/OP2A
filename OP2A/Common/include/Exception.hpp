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

class Common_API ExceptionManager : public Common::NonCopyable<ExceptionManager>
{
public:
	ExceptionManager();

	static ExceptionManager& getInstance ();

	bool	ExceptionOutputs;
	bool	ExceptionDumps;
	bool	ExceptionAborts;
};


class Common_API Exception : public std::exception
{
public:
	  virtual ~Exception () throw ();

	  std::string full_description ()	const throw ();
	  const std::string& str () 		const throw ();

	  const char* what () const throw ();
	  void append (const std::string& add) throw ();

	  virtual std::string getClassName () throw () { return m_class_name; }

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
