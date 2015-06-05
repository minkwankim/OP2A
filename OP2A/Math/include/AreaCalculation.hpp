/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 5, 2015
 *      			Author: Minkwan Kim
 *
 * AreaCalculation.hpp
 * 			-  
 *  
 */
#ifndef AREACALCULATION_HPP_
#define AREACALCULATION_HPP_

#include "Math/include/OP2A_Math.hpp"

namespace OP2A{
namespace Math{



// AREA/LENGTH FUNCTIONS
double length(std::vector<double> n1, std::vector<double> n2);
double CalAreaTriangle(std::vector<double> n1, std::vector<double> n2,  std::vector<double> n3);
double CalAreaQuadrilateral(std::vector<double> n1, std::vector<double> n2,  std::vector<double> n3, std::vector<double> n4);
double CalVolumeTetrahedron(std::vector<double> n1, std::vector<double> n2,  std::vector<double> n3, std::vector<double> n4);
double CalVolumePyramid(std::vector<double> n1, std::vector<double> n2,  std::vector<double> n3, std::vector<double> n4,  std::vector<double> n5);
double CalVolumeWedge(std::vector<double> n1, std::vector<double> n2,  std::vector<double> n3, std::vector<double> n4, std::vector<double> n5, std::vector<double> n6);
double CalVolumeHexahedron(std::vector<double> n1, std::vector<double> n2,  std::vector<double> n3, std::vector<double> n4, std::vector<double> n5, std::vector<double> n6, std::vector<double> n7, std::vector<double> n8);




}
}




#endif /* AREACALCULATION_HPP_ */
