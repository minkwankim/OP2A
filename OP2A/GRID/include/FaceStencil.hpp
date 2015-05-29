/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 29, 2015
 *      			Author: Minkwan Kim
 *
 * FaceStencil.hpp
 * 			-  
 *  
 */
#ifndef FACESTENCIL_HPP_
#define FACESTENCIL_HPP_


#include <iostream>
#include <vector>

using namespace std;


namespace OP2A{
namespace GRID{

class Cell;

enum StencilLabel
{
	cl			= 0,
	cll			= 1,
	clll		= 2,
	clupper		= 3,
	cllower		= 4,

	cr			= 0,
	crr			= 1,
	crrr		= 2,
	crupper		= 3,
	crlower		= 4,
};




/*
 * Class for face stencil information
 * @author Minkwan Kim
 * @version 1.0 27/5/2015
 */
class FaceStencil
{
public:
	vector<Cell*>	cl;
	vector<Cell*>	cr;

	FaceStencil();
	explicit FaceStencil(bool extendedStencil);

	void apply_extended_sencil();
	bool is_extended();

	~FaceStencil();

private:
	bool extended;
};



}
}


#endif /* FACESTENCIL_HPP_ */
