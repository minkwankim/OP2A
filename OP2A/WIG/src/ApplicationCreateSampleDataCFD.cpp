/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 11, 2015
 *      			Author: Minkwan Kim
 *
 * ApplicationCreateSampleDataCFD.cpp
 * 			-  
 *  
 */



#include "Common/include/CodeLocation.hpp"
#include "Common/include/Exception_NPExceed.hpp"


#include "../include/OP2A_Application.hpp"






void ApplicationOP2A::create_sampleDataCFD()
{
	int numData;
	numData	= species_set.NS;
	numData	+= grid.ND;
	numData	+= 1;
	if (problem_setup.NER != 0) numData	+= 1;
	if (problem_setup.NEV != 0) numData	+= 1;
	if (problem_setup.NEE != 0) numData	+= 1;

	Common::Map1D<std::string, int>	QnameMap (numData);
	Common::Map1D<std::string, int>	VnameMap (numData);
	Common::Map1D<std::string, int>	WnameMap (numData);



	for (int s = 0; s <= species_set.NS-1; s++)
	{
		QnameMap.insert("rho<sub>" + species_set.species[s].name + "</sub> [kg/m<sup>-3</sup>]", s);
		VnameMap.insert("rho<sub>" + species_set.species[s].name + "</sub> [kg/m<sup>-3</sup>]", s);
		WnameMap.insert("rho<sub>" + species_set.species[s].name + "</sub> [kg/m<sup>-3</sup>]", s);
	}


	QnameMap.insert("rho u", species_set.NS);
	VnameMap.insert("u [m/s]", species_set.NS);
	WnameMap.insert("u [m/s]", species_set.NS);

	QnameMap.insert("rho v", species_set.NS+1);
	VnameMap.insert("v [m/s]", species_set.NS+1);
	WnameMap.insert("v [m/s]", species_set.NS+1);


	if (grid.ND == 3)
	{
		QnameMap.insert("rho w", species_set.NS+2);
		VnameMap.insert("w [m/s]", species_set.NS+2);
		WnameMap.insert("w [m/s]", species_set.NS+2);
	}

	QnameMap.insert("E <sub>total</sub> [J]", species_set.NS+grid.ND);
	VnameMap.insert("T <sub>tra</sub> [K]", species_set.NS+grid.ND);
	WnameMap.insert("P [N/m<sup>2</sup>]", species_set.NS+grid.ND);

	int ne	= 0;
	if (problem_setup.NER != 0)
	{
		ne	= ne +1;
		QnameMap.insert("E <sub>rot</sub> [J]", species_set.NS+grid.ND+ne);
		VnameMap.insert("T <sub>rot</sub> [K]", species_set.NS+grid.ND+ne);
		WnameMap.insert("T <sub>rot</sub> [K]", species_set.NS+grid.ND+ne);
	}

	if (problem_setup.NEV != 0)
	{
		ne	= ne + 1;
		QnameMap.insert("E <sub>vib</sub> [J]", species_set.NS+grid.ND+ne);
		VnameMap.insert("T <sub>vib</sub> [K]", species_set.NS+grid.ND+ne);
		WnameMap.insert("T <sub>vib</sub> [K]", species_set.NS+grid.ND+ne);
	}

	if (problem_setup.NEE != 0)
	{
		ne	= ne + 1;
		QnameMap.insert("E <sub>e</sub> [J]", species_set.NS+grid.ND+ne);
		VnameMap.insert("T <sub>e</sub> [K]", species_set.NS+grid.ND+ne);
		WnameMap.insert("T <sub>e</sub> [K]", species_set.NS+grid.ND+ne);
	}


	data_CFD_Q.resize(numData, QnameMap);	data_CFD_Q.asgName("Conservative variables");
	data_CFD_V.resize(numData, VnameMap);	data_CFD_V.asgName("Primitive variables");
	data_CFD_W.resize(numData, WnameMap);	data_CFD_W.asgName("W");


	int nCv	= 1;
	if (problem_setup.NER != 0)	nCv++;

	Common::Map1D<std::string, int>	MixtureNameMap (4+nCv);
	MixtureNameMap.insert("rho_mix", 0);
	MixtureNameMap.insert("R_mix", 1);
	MixtureNameMap.insert("M_mix", 2);
	MixtureNameMap.insert("gamma", 3);
	if (problem_setup.NER != 0)
	{
		MixtureNameMap.insert("Cv_tra_mix", 4);
		MixtureNameMap.insert("Cv_rot_mix", 5);
	}
	else
	{
		MixtureNameMap.insert("Cv_tr_mix", 4);
	}
	data_CFD_mixture.resize(4+nCv, MixtureNameMap); data_CFD_Q.asgName("Mixture variables");


	Common::Map1D<std::string, int>	YsNameMap (species_set.NS);
	Common::Map1D<std::string, int>	XsNameMap (species_set.NS);

	for (int s = 0; s <= species_set.NS-1; s++)
	{
		YsNameMap.insert("Y<sub>" + species_set.species[s].name + "</sub>", s);
		XsNameMap.insert("X<sub>" + species_set.species[s].name + "</sub>", s);
	}
	data_CFD_Xs.resize(species_set.NS, XsNameMap); data_CFD_Xs.asgName("Species mole fraction");
	data_CFD_Ys.resize(species_set.NS, YsNameMap); data_CFD_Ys.asgName("Species mass fraction");







}
