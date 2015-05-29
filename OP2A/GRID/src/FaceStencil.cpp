/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 29, 2015
 *      			Author: Minkwan Kim
 *
 * FaceStencil.cpp
 * 			-  
 *  
 */



#include "../include/FaceStencil.hpp"

namespace OP2A{
namespace GRID{


FaceStencil::FaceStencil(bool extendedStencil)
{
	if (extendedStencil == true)
	{
		cl.resize(5);
		cr.resize(5);
	}
	else
	{
		cl.resize(2);
		cr.resize(2);
	}
}


FaceStencil::~FaceStencil()
{

}



}
}
