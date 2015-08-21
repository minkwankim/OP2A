/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 22, 2015
 *      			Author: Minkwan Kim
 *
 * BoundaryConditions.hpp
 * 			-  
 *  
 */
#ifndef BOUNDARYCONDITIONS_HPP_
#define BOUNDARYCONDITIONS_HPP_

#include "GRID/include/GridConstant.hpp"
#include "DATA/include/DataStorage.hpp"
#include "DATA/include/DataStorage2D.hpp"

#include "CFD/include/VariableChange.hpp"
#include "CFD/include/VariableSet.hpp"
#include "CFD/include/VariableConstants.hpp"


namespace OP2A{
namespace CFD{

enum CFDBCTypes
{
	Interior		= 0,
	WallType		= 1,
	InletType		= 2,
	FreestreamType 	= 3,
	ExitType		= 4,
	OtherType		= 5
};


double CFD_CalculateWallTemperature(double Tcl, double kappa_cl, double dn, bool use_emissivity, double Tref, vector<double> emissivity_below, vector<double> emissivity_above, int MaxIter, double epsilon);



class CFD_API BCInviscid: public Common::NonInstantiable<BCInviscid>
{
public:
	static CFDBCTypes	BCTypeInCFD(const int FaceBCType);

	static void  wallTypeBC(Data::DataStorage& Qcl, vector< vector<double> >& face_normal_vectors, Data::DataStorage& Qcr, int ND);
	static void inletTypeBC(Data::DataStorage& Qcl, Data::DataStorage& Qcr, int ND, Data::DataStorage& Qinlet);
	static void  exitTypeBC(Data::DataStorage& Qcl, Data::DataStorage& Qcr, int ND);
};

class CFD_API BCInviscidImplicit: public Common::NonInstantiable<BCInviscid>
{
public:
	static void  wallTypeBC(Data::DataStorage2D& J_plus, Data::DataStorage2D& J_minus, vector< vector<double> >& face_normal_vector, int NS, int ND, int NE);
	static void  axisTypeBC(Data::DataStorage2D& J_plus, Data::DataStorage2D& J_minus, vector< vector<double> >& face_normal_vector, int NS, int ND, int NE);

	static void inletTypeBC(Data::DataStorage2D& J_plus, Data::DataStorage2D& J_minus, vector< vector<double> >& face_normal_vector, int NS, int ND, int NE);
	static void  exitTypeBC(Data::DataStorage2D& J_plus, Data::DataStorage2D& J_minus, vector< vector<double> >& face_normal_vector, int NS, int ND, int NE);
};



class CFD_API BCViscous: public Common::NonInstantiable<BCViscous>
{
public:
	static CFDBCTypes	BCTypeInCFD(const int FaceBCType);

	static void  wallTypeBC(Data::DataStorageVector<Data::DataStorage>& CellData1D_cl,
							Data::DataStorageVector<Data::DataStorage>& CellData1D_cr,
							CHEM::SpeciesSet& species_set,
							Data::DataStorage& BCValuesWall,
							vector<double>& xcl,
							vector<double>& xf,
							int  ND,
							int  CFD_variabletype,
							int	 CFD_NT,
							bool adiabaticWall,
							bool catalyticWall,
							bool nonSlipWall,
							bool radiativeWall,
							double kappa_cl,
							bool use_emissivity,
							double Tref,
							vector<double> emissivity_below,
							vector<double> emissivity_above);

	static void OtherTypeBC(Data::DataStorageVector<Data::DataStorage>& CellData1D_cl, Data::DataStorageVector<Data::DataStorage>& CellData1D_cr, int ND);
};


class CFD_API BCViscousImplicit: public Common::NonInstantiable<BCViscousImplicit>
{
public:
	static void  wallTypeBC(Data::DataStorage& Qcl, vector<double>& face_normal_vector, Data::DataStorage& Qcr, int ND);
	static void inletTypeBC(Data::DataStorage& Qcl, Data::DataStorage& Qcr, int ND, Data::DataStorage& Qinlet);
	static void  exitTypeBC(Data::DataStorage& Qcl, Data::DataStorage& Qcr, int ND);
};



}
}

#endif /* BOUNDARYCONDITIONS_HPP_ */
