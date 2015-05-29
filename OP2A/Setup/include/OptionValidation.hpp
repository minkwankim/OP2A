/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 26, 2015
 *      			Author: Minkwan Kim
 *
 * OptionValidation.hpp
 * 			-  
 *  
 */
#ifndef OPTIONVALIDATION_HPP_
#define OPTIONVALIDATION_HPP_

#include "Common/include/OP2A.hpp"
#include "Setup/include/SetupAPI.hpp"


namespace OP2A{
namespace Setup{

class Option;

/*
 * Class to validate a given option based on some restriction defined in the derived classes
 * @author	Minkwan Kim
 * @version 1.0 25/05/2015
 */
class Setup_API OptionValidation
{
public:
	OptionValidation (Option * opt);
	virtual ~OptionValidation();

	virtual bool isValid()	= 0;
	virtual std::string getReason() = 0;


protected:
	Option * m_option;	// ptr to the concrete option to validate

};


}
}





#endif /* OPTIONVALIDATION_HPP_ */
