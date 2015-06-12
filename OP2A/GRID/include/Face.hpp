/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 1, 2015
 *      			Author: Minkwan Kim
 *
 * Face.hpp
 * 			-  
 *  
 */
#ifndef FACE_HPP_
#define FACE_HPP_



#include <vector>
using namespace std;

#include "FaceGeometry.hpp"

#include "DATA/include/DataStorage.hpp"
#include "DATA/include/DataStorage2D.hpp"
#include "DATA/include/DataStorageVector.hpp"

namespace OP2A{
namespace GRID{


/*
 * Class for Face
 * @author Minkwan Kim
 * @version 1.0 27/5/2015
 */

class Face
{
public:
	FaceGeometry									geo;
	Data::DataStorageVector<Data::DataStorage>		data1D;
	Data::DataStorageVector<Data::DataStorage2D>	data2D;


	Face():allocated_geo(false), allocated_data1D(false), allocated_data2D(false)
	{

	}

	explicit	Face(const int ND, const FaceType f_type, bool extendedStencil,	Data::DataStorageVector<Data::DataStorage>& data_template);
	explicit	Face(FaceGeometry& geo_template, Data::DataStorageVector<Data::DataStorage>& data_template);

	explicit	Face(const int ND, const FaceType f_type, bool extendedStencil,	Data::DataStorageVector<Data::DataStorage>& data1D_template, Data::DataStorageVector<Data::DataStorage2D>& data2D_template);
	explicit	Face(FaceGeometry& geo_template, Data::DataStorageVector<Data::DataStorage>& data1D_template, Data::DataStorageVector<Data::DataStorage2D>& data2D_template);

	~Face();

public:
	void allocate_geo_type(const int ND, const FaceType f_type, bool extendedStencil);
	void allocate_geo_type(FaceGeometry& geo_template);
	void allocate_data1D_type(Data::DataStorageVector<Data::DataStorage>& data1D_template);
	void allocate_data2D_type(Data::DataStorageVector<Data::DataStorage2D>& data2D_template);

	bool is_allocated_geo();
	bool is_allocated_data1D();
	bool is_allocated_data2D();


private:
	bool allocated_geo;
	bool allocated_data1D;
	bool allocated_data2D;
};








}
}
#endif /* FACE_HPP_ */
