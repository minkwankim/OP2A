/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 26, 2015
 *      			Author: Minkwan Kim
 *
 * ProblemBasic.hpp
 * 			-  
 *  
 */
#ifndef PROBLEMBASIC_HPP_
#define PROBLEMBASIC_HPP_

#include "Setup/include/SetupObject.hpp"


using namespace std;

namespace OP2A{
namespace Problem{


class ProblemBasic: public OP2A::Setup::SetupObject
{
public:
	string	title;
	OPuint	module;




};




}
}




#endif /* PROBLEMBASIC_HPP_ */
