/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 21, 2015
 *      			Author: Minkwan Kim
 *
 * ProblemVars.hpp
 * 			-  
 *  
 */
#ifndef PROBLEMVARS_HPP_
#define PROBLEMVARS_HPP_


#include	"Common/include/OP2A.hpp"
#include	"Common/include/NonCopyable.hpp"

#include	"Problem/include/ProblemAPI.hpp"

namespace OP2A{
namespace Problem{




class Problem_API	ProblemVars: public Common::NonCopyable<ProblemVars>
{
public:
	ProblemVars();


};

}
}


#endif /* PROBLEMVARS_HPP_ */
