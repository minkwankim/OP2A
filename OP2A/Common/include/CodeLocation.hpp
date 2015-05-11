/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 11, 2015
 *      			Author: Minkwan Kim
 *
 * CodeLocation.hpp
 * 			-  
 *  
 */
#ifndef CODELOCATION_HPP_
#define CODELOCATION_HPP_


#include "Common/include/Common.hpp"

namespace OP2A {
namespace Common {

/*
 * This class stores the information about a location in the source code
 */

class Common_API Code_location
{
	public:
	  explicit Code_location (const char *file, int line, const char *function);

	  std::string str () const;

	private:
	  std::string	m_file;
	  std::string	m_function;
	  int			m_line;
};
//////////////////////////////////////////////////////////////////////////////


#define FromHere() OP2A::Common::Code_location( __FILE__ , __LINE__ , __FUNCTION__ )

}
}



#endif /* CODELOCATION_HPP_ */
