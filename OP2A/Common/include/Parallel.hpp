/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 14, 2015
 *      			Author: Minkwan Kim
 *
 * Parallel.hpp
 * 			-  
 *  
 */
#ifndef PARALLEL_HPP_
#define PARALLEL_HPP_

#include <time.h>

#include "Common/include/Standard_headers.hpp"
#include "Common/include/CommonAPI.hpp"


//////////////////////////////////////////////////////////////////////////////
namespace OP2A {
namespace Common{

/*
 * Class for 2D array (Matrix, I x J)
 * 	- [Size_]: size of Coloumn (I)
 * 	- [ElementSize_]: Size of row (J)
 * 	- [Data_]: Data of matrix
 *
 * 	Initially written by:	Minkwan Kim
 * 	Last modified on	:	May/14/2015
 * 	Last modified by	:	Minkwan Kim
 */
class Parallel_info
{
public:
	time_t	t;
	struct	tm * now;

	/*
	OPdouble 	t0;
	OPdouble 	StartTime;
	OPdouble	StopTime;
	OPdouble 	UsedTime;

	OPint		rank;
	OPint		size;
*/
	Parallel_info();
	~Parallel_info();

};



}
}


#endif /* PARALLEL_HPP_ */
