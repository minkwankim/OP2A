/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 14, 2015
 *      			Author: Minkwan Kim
 *
 * Exception_Universial.hpp
 * 			-  
 *  
 */
#ifndef EXCEPTION_UNIVERSIAL_HPP_
#define EXCEPTION_UNIVERSIAL_HPP_

#include "Common/include/Exception.hpp"

namespace OP2A {
namespace Common {


class Common_API Exception_Universial : public Common::Exception
{
public:
	Exception_Universial (const Common::Code_location& where, const std::string& what, const std::string& reason) : Common::Exception(where, what, reason) {}

	/// Copy constructor
	Exception_Universial(const Exception_Universial& e) throw  () : Exception(e) {}
};


}
}


#endif /* EXCEPTION_UNIVERSIAL_HPP_ */
