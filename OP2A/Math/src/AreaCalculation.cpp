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

#include <limits>
#include "Math/include/AreaCalculation.hpp"
#include "Math/include/OP2A_Vector.hpp"

#include "Common/include/Exception_DimensionMatch.hpp"
#include "Common/include/Exception_NaNValue.hpp"
#include "Common/include/Exception_NegativeValue.hpp"
#include "Common/include/Exception_InfiniteValue.hpp"




using namespace std;
namespace OP2A{
namespace Math{




double length(std::vector<double> n1, std::vector<double> n2)
{
	int N = n1.size();
	if (N != n2.size())	throw Common::ExceptionDimensionMatch (FromHere(), "Dimenction of X and Y values do not match");

	double len	= 0.0;

	for (int i = 0; i <= N-1; i++)	len	+= pow(n1[i] - n2[i], 2.0);

	return (sqrt(len));
}



double CalAreaTriangle(std::vector<double> n1, std::vector<double> n2,  std::vector<double> n3)
{
	double a, b, c, s;
	double area;

	a = length(n1, n2);
	b = length(n2, n3);
	c = length(n3, n1);

	s = 0.5 * (a + b + c);

	area	= s * (s-a) * (s-b) * (s-c);
	if (area < 0.0)	throw Common::ExceptionNegativeValue (FromHere(), "Negative value in the calculation of triangle ares");

	return (sqrt(area));
}



double CalAreaQuadrilateral(std::vector<double> n1, std::vector<double> n2,  std::vector<double> n3, std::vector<double> n4)
{
	double s1, s2;

	s1 = CalAreaTriangle(n1, n2, n3);
	s2 = CalAreaTriangle(n2, n3, n4);

	return (s1 + s2);
}


double CalVolumeTetrahedron(std::vector<double> n1, std::vector<double> n2,  std::vector<double> n3, std::vector<double> n4)
{
	double vol;
	VECTOR A(n4, n1);
	VECTOR B(n4, n2);
	VECTOR C(n4, n3);

	VECTOR TEMP1	= VectorCrossProduct(B, C);
	vol		= VectorDotProduct(A, TEMP1);
	vol		= fabs(vol) / 6.0;

	if (vol	!= vol) throw Common::ExceptionNaNValue (FromHere(), "Nan value for Tetrahedron volumen ");
	if (vol == numeric_limits<double>::infinity()) throw Common::ExceptionInfiniteValue (FromHere(), "Inifnite Value for Tetrahedron volumen");

	return (vol);
}


double CalVolumePyramid(std::vector<double> n1, std::vector<double> n2,  std::vector<double> n3, std::vector<double> n4,  std::vector<double> n5)
{
	int i;
	double volume;
	std::vector<double> n6(3, 0.0);


	for (i = 0; i <= 2; i++)	n6[i] = 0.25 * (n1[i] + n2[i] + n3[i] + n4[i]);


	volume = CalVolumeTetrahedron(n1, n2, n6, n5);
	volume += CalVolumeTetrahedron(n2, n3, n6, n5);
	volume += CalVolumeTetrahedron(n3, n4, n6, n5);
	volume += CalVolumeTetrahedron(n4, n1, n6, n5);

	if (volume	!= volume) throw Common::ExceptionNaNValue (FromHere(), "Nan value for Tetrahedron volumen ");
	if (volume == numeric_limits<double>::infinity()) throw Common::ExceptionInfiniteValue (FromHere(), "Inifnite Value for Tetrahedron volumen");

	return (volume);
}


double CalVolumeWedge(std::vector<double> n1, std::vector<double> n2,  std::vector<double> n3, std::vector<double> n4, std::vector<double> n5, std::vector<double> n6)
{
	int i;
	std::vector<double> n7(3, 0.0);

	double volume;

	for (i = 0; i <= 2; i++)	n7[i] = (n1[i] + n2[i] + n3[i] + n4[i] + n5[i] + n6[i]) / 6.0;

	volume = CalVolumePyramid(n1, n4, n5, n2, n7);
	volume += CalVolumePyramid(n2, n3, n6, n3, n7);
	volume += CalVolumePyramid(n1, n3, n6, n4, n7);
	volume += CalVolumeTetrahedron(n1, n2, n3, n7);
	volume += CalVolumeTetrahedron(n4, n6, n5, n7);

	if (volume	!= volume) throw Common::ExceptionNaNValue (FromHere(), "Nan value for Tetrahedron volumen ");
	if (volume == numeric_limits<double>::infinity()) throw Common::ExceptionInfiniteValue (FromHere(), "Inifnite Value for Tetrahedron volumen");

	return (volume);
}



double CalVolumeHexahedron(std::vector<double> n1, std::vector<double> n2,  std::vector<double> n3, std::vector<double> n4, std::vector<double> n5, std::vector<double> n6, std::vector<double> n7, std::vector<double> n8)
{
	int i;
	std::vector<double> n9(3, 0.0);
	double volume;

	for (i = 0; i <= 2; i++)	n9[i] = 0.125 * (n1[i] + n2[i] + n3[i] + n4[i] + n5[i] + n6[i] + n7[i] + n8[i]);


	volume = CalVolumePyramid(n1, n2, n3, n4, n9);
	volume += CalVolumePyramid(n5, n8, n7, n6, n9);
	volume += CalVolumePyramid(n1, n5, n6, n2, n9);
	volume += CalVolumePyramid(n3, n7, n8, n4, n9);
	volume += CalVolumePyramid(n2, n6, n7, n3, n9);
	volume += CalVolumePyramid(n1, n4, n8, n5, n9);

	if (volume	!= volume) throw Common::ExceptionNaNValue (FromHere(), "Nan value for Tetrahedron volumen ");
	if (volume == numeric_limits<double>::infinity()) throw Common::ExceptionInfiniteValue (FromHere(), "Inifnite Value for Tetrahedron volumen");

	return (volume);
}





}
}
