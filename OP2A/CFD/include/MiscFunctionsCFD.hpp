/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 24, 2015
 *      			Author: Minkwan Kim
 *
 * MiscFunctionsCFD.cpp
 * 			-  
 *  
 */
#ifndef MISCFUNCTIONSCFD_CPP_
#define MISCFUNCTIONSCFD_CPP_


namespace OP2A{
namespace CFD{

	double CalculateCFLNumber(const int iter, const int iterBefore1, const double CFLStart, const double CFLMax);
	double Calculatedt(const double dX, const double U, const double Maxdt, const double CFL);




}
}
#endif /* MISCFUNCTIONSCFD_CPP_ */
