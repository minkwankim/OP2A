/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Sep 1, 2015
 *      			Author: Minkwan Kim
 *
 * GradientCalculation.hpp
 * 			-  
 *  
 */
#ifndef GRADIENTCALCULATION_HPP_
#define GRADIENTCALCULATION_HPP_


#include <limits>
#include "Common/include/Exception_NaNValue.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"


#include "CFD/include/CFDCommon.hpp"
#include "DATA/include/DataStorage2D.hpp"
#include "DATA/include/DataStorageVector.hpp"

#include "GRID/include/Grid.hpp"



namespace OP2A{
namespace CFD{


class CFD_API GradientCalculation : public Common::NonInstantiable<GradientCalculation>
{
public:
	static double GradientCalculation_FaceAveraging(double phi_cl, double phi_cr);
	static double GradientCalculation_IDW(vector<double>& x_f, vector< vector<double> >& x_i, vector <double>& phi_i, int N, int ND);
	//static double GradientCalculation_LSQR(vector<double>& x_f, vector< vector<double> >& x_i, vector <double>& phi_i, int N, int ND);

	static void GradientCalculation_GG(GRID::Grid &grid, int i_VAR1, int i_VAR2, int method, vector<vector <double> >& grad_phi);
	static void GradientCalculation_LSQR(GRID::Grid &grid, int i_VAR1, int i_VAR2, int method, vector<vector <double> >& grad_phi);


};





}
}


#endif /* GRADIENTCALCULATION_HPP_ */
