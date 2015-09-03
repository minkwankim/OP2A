/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Sep 1, 2015
 *      			Author: Minkwan Kim
 *
 * GradientCalculation_FaceAveraging.cpp
 * 			-  
 *  
 */



#include <omp.h>
#include <time.h>

#include "CFD/include/GradientCalculation.hpp"

#include "Math/include/OP2A_Matrix.hpp"
#include "Math/include/MathMisc.hpp"
#include "Common/include/Time_Info.hpp"
#include "Common/include/MultiDimension.hpp"

namespace OP2A{
namespace CFD{

double GradientCalculation::GradientCalculation_FaceAveraging(double phi_cl, double phi_cr)
{
	double phi_f;
	phi_f	= (phi_cl + phi_cr) / 2.0;
	return (phi_f);
}





}
}
