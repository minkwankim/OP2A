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
#include "Common/include/StringOps.hpp"

#include "../include/OP2A_Application.hpp"
#include "../include/ApplicationConstants.hpp"






void ApplicationOP2A::create_sampleDataCFD()
{
	int numData;
	numData	= species_set.NS;
	numData	+= grid.ND;
	numData	+= 1;
	if (problem_setup.NER != 0) numData	+= 1;
	if (problem_setup.NEV != 0) numData	+= 1;
	if (problem_setup.NEE != 0) numData	+= 1;

	vector<string>	Q(numData);
	vector<string>	V(numData);
	vector<string>	F(numData);


	/*
	 * Assign Sample data format for
	 * 	 - Q
	 * 	 - V
	 * 	 - W
	 * 	 - Flux (viscous/inviscid)
	 */
	Common::Map1D<std::string, int>	QnameMap (numData);
	Common::Map1D<std::string, int>	VnameMap (numData);
	Common::Map1D<std::string, int>	WnameMap (numData);
	Common::Map1D<std::string, int>	FluxInvnameMap (numData);
	Common::Map1D<std::string, int>	FluxVisnameMap (numData);

	// Species
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Q[s]	= string(VAR_RHO) + species_set.species[s].name;
		V[s]	= string(VAR_RHO) + species_set.species[s].name;
		F[s]	= string(VAR_FLUX) + string(VAR_RHO) + species_set.species[s].name;

		QnameMap.insert(string(VAR_RHO_PRE) + species_set.species[s].name + string(VAR_RHO_POST), s);
		VnameMap.insert(string(VAR_RHO_PRE) + species_set.species[s].name + string(VAR_RHO_POST), s);
		WnameMap.insert(string(VAR_RHO_PRE) + species_set.species[s].name + string(VAR_RHO_POST), s);

		FluxInvnameMap.insert(string(VAR_F_INV) + string(VAR_SPECIES) + Common::StringOps::to_str<int>(s), s);
		FluxVisnameMap.insert(string(VAR_F_VIS) + string(VAR_SPECIES) + Common::StringOps::to_str<int>(s), s);
	}

	// Momentum
	//  - X
	Q[species_set.NS]	= VAR_RHO_U;
	V[species_set.NS]	= VAR_U_WO_UNIT;
	F[species_set.NS]	= string(VAR_FLUX) + string(VAR_RHO_U);
	QnameMap.insert(VAR_RHO_U, species_set.NS);
	VnameMap.insert(VAR_U, species_set.NS);
	WnameMap.insert(VAR_U, species_set.NS);

	FluxInvnameMap.insert(string(VAR_F_INV)+string(VAR_RHO_U), species_set.NS);
	FluxVisnameMap.insert(string(VAR_F_VIS)+string(VAR_RHO_U), species_set.NS);

	//  - Y
	Q[species_set.NS+1]	= VAR_RHO_V;
	V[species_set.NS+1]	= VAR_V_WO_UNIT;
	F[species_set.NS+1]	= string(VAR_FLUX) + string(VAR_RHO_V);
	QnameMap.insert(VAR_RHO_V, species_set.NS+1);
	VnameMap.insert(VAR_V, species_set.NS+1);
	WnameMap.insert(VAR_V, species_set.NS+1);

	FluxInvnameMap.insert(string(VAR_F_INV)+string(VAR_RHO_V), species_set.NS+1);
	FluxVisnameMap.insert(string(VAR_F_VIS)+string(VAR_RHO_V), species_set.NS+1);

	//  - Z
	if (grid.ND == 3)
	{
		Q[species_set.NS+2]	= VAR_RHO_W;
		V[species_set.NS+2]	= VAR_W_WO_UNIT;
		F[species_set.NS+2]	= string(VAR_FLUX) + string(VAR_RHO_W);

		QnameMap.insert(VAR_RHO_W, species_set.NS+2);
		VnameMap.insert(VAR_W, species_set.NS+2);
		WnameMap.insert(VAR_W, species_set.NS+2);

		FluxInvnameMap.insert(string(VAR_F_INV)+string(VAR_RHO_W), species_set.NS+2);
		FluxVisnameMap.insert(string(VAR_F_VIS)+string(VAR_RHO_W), species_set.NS+2);
	}


	// Total Energy
	Q[species_set.NS+grid.ND]	= VAR_E_WO_UNIT;
	V[species_set.NS+grid.ND]	= VAR_T_WO_UNIT;
	F[species_set.NS+grid.ND]	= string(VAR_FLUX) + string(VAR_E_WO_UNIT);;
	QnameMap.insert(VAR_E, species_set.NS+grid.ND);
	VnameMap.insert(VAR_T, species_set.NS+grid.ND);
	WnameMap.insert(VAR_P, species_set.NS+grid.ND);

	int ne	= 0;
	FluxInvnameMap.insert(string(VAR_F_INV)+string(VAR_E_WO_UNIT), species_set.NS+grid.ND+ne);
	FluxVisnameMap.insert(string(VAR_F_VIS)+string(VAR_E_WO_UNIT), species_set.NS+grid.ND+ne);
	ne++;

	// Rotational energy
	if (problem_setup.NER != 0)
	{
		Q[species_set.NS+grid.ND+ne]	= VAR_E_ROT_WO_UNIT;
		V[species_set.NS+grid.ND+ne]	= VAR_T_ROT_WO_UNIT;
		F[species_set.NS+grid.ND+ne]	= string(VAR_FLUX) + string(VAR_E_ROT_WO_UNIT);
		QnameMap.insert(VAR_E_ROT, species_set.NS+grid.ND+ne);
		VnameMap.insert(VAR_T_ROT, species_set.NS+grid.ND+ne);
		WnameMap.insert(VAR_T_ROT, species_set.NS+grid.ND+ne);

		FluxInvnameMap.insert(string(VAR_F_INV) + string(VAR_E_ROT_WO_UNIT), species_set.NS+grid.ND+ne);
		FluxVisnameMap.insert(string(VAR_F_INV) + string(VAR_E_ROT_WO_UNIT), species_set.NS+grid.ND+ne);
		ne++;
	}

	if (problem_setup.NEV != 0)
	{
		Q[species_set.NS+grid.ND+ne]	= VAR_E_VIB_WO_UNIT;
		V[species_set.NS+grid.ND+ne]	= VAR_T_VIB_WO_UNIT;
		F[species_set.NS+grid.ND+ne]	= string(VAR_FLUX) + string(VAR_E_VIB_WO_UNIT);
		QnameMap.insert(VAR_E_VIB, species_set.NS+grid.ND+ne);
		VnameMap.insert(VAR_T_VIB, species_set.NS+grid.ND+ne);
		WnameMap.insert(VAR_T_VIB, species_set.NS+grid.ND+ne);

		FluxInvnameMap.insert(string(VAR_F_INV) + string(VAR_E_VIB_WO_UNIT), species_set.NS+grid.ND+ne);
		FluxVisnameMap.insert(string(VAR_F_VIS) + string(VAR_E_VIB_WO_UNIT), species_set.NS+grid.ND+ne);

		ne++;
	}

	if (problem_setup.NEE != 0)
	{
		Q[species_set.NS+grid.ND+ne]	= VAR_E_E_WO_UNIT;
		V[species_set.NS+grid.ND+ne]	= VAR_T_E_WO_UNIT;
		F[species_set.NS+grid.ND+ne]	= string(VAR_FLUX) + string(VAR_E_E_WO_UNIT);
		QnameMap.insert(VAR_E_E, species_set.NS+grid.ND+ne);
		VnameMap.insert(VAR_T_E, species_set.NS+grid.ND+ne);
		WnameMap.insert(VAR_T_E, species_set.NS+grid.ND+ne);

		FluxInvnameMap.insert(string(VAR_F_INV) + string(VAR_E_E_WO_UNIT), species_set.NS+grid.ND+ne);
		FluxVisnameMap.insert(string(VAR_F_VIS) + string(VAR_E_E_WO_UNIT), species_set.NS+grid.ND+ne);
		ne++;
	}


	data_CFD_Q.resize(numData, QnameMap);	data_CFD_Q.asgName(VAR_VECTOR_Q);
	data_CFD_V.resize(numData, VnameMap);	data_CFD_V.asgName(VAR_VECTOR_V);
	data_CFD_W.resize(numData, WnameMap);	data_CFD_W.asgName(VAR_VECTOR_W);

	data_CFD_Flux_inviscid.resize(numData, FluxInvnameMap);	data_CFD_Flux_inviscid.asgName(VAR_VECTOR_F_INV);
	data_CFD_Flux_viscous.resize(numData, FluxVisnameMap);	data_CFD_Flux_viscous.asgName(VAR_VECTOR_F_VIS);



	/*
	 * Assign Sample data format for
	 * 	 - Mixture data
	 */

	int nCv	= 1;
	if (problem_setup.NER != 0)	nCv++;

	int nMixture;
	if (problem_setup.is_viscous == true)	nMixture	= 4+nCv+1+ne;
	else									nMixture	= 4+nCv;

	Common::Map1D<std::string, int>	MixtureNameMap (nMixture);

	MixtureNameMap.insert(string(VAR_RHO) + string(VAR_MIX), 0);
	MixtureNameMap.insert(string(VAR_R)	+ string(VAR_MIX_SUB), 1);
	MixtureNameMap.insert(string(VAR_M)	+ string(VAR_MIX_SUB), 2);
	MixtureNameMap.insert(string(VAR_GAMMA)	+ string(VAR_MIX_SUB), 3);

	if (problem_setup.NER != 0)
	{
		MixtureNameMap.insert(string(VAR_CV_TRA) + string(VAR_MIX_SUB), 4);
		MixtureNameMap.insert(string(VAR_CV_ROT) + string(VAR_MIX_SUB), 5);
	}
	else
	{
		MixtureNameMap.insert(string(VAR_CV_TR) + string(VAR_MIX_SUB), 4);
	}

	if (problem_setup.is_viscous == true)
	{
		MixtureNameMap.insert(string(VAR_VISCOSITY_COEFF) + string(VAR_MIX_SUB), 4+nCv);
		ne	= 0;

		MixtureNameMap.insert(string(VAR_THERMAL_COND) + string(VAR_TRA_SUB), 4+nCv+ne);
		ne++;

		if (problem_setup.NER != 0)
		{
			MixtureNameMap.insert(string(VAR_THERMAL_COND) + string(VAR_ROT_SUB), 4+nCv+ne);
			ne++;
		}

		if (problem_setup.NEV != 0)
		{
			MixtureNameMap.insert(string(VAR_THERMAL_COND) + string(VAR_VIB_SUB), 4+nCv+ne);
			ne++;
		}

		if (problem_setup.NEE != 0)
		{
			MixtureNameMap.insert(string(VAR_THERMAL_COND) + string(VAR_E_SUB), 4+nCv+ne);
			ne++;
		}
	}

	data_CFD_mixture.resize(nMixture, MixtureNameMap); data_CFD_Q.asgName(VAR_VECTOR_MIX);



	/*
	 * Assign Sample data format for
	 * 	 - Xs
	 * 	 - Ys
	 */
	Common::Map1D<std::string, int>	YsNameMap (species_set.NS);
	Common::Map1D<std::string, int>	XsNameMap (species_set.NS);

	for (int s = 0; s <= species_set.NS-1; s++)
	{
		YsNameMap.insert(string(VAR_YS) + species_set.species[s].name, s);
		XsNameMap.insert(string(VAR_XS) + species_set.species[s].name, s);
	}

	data_CFD_Xs.resize(species_set.NS, XsNameMap); data_CFD_Xs.asgName(VAR_VECTOR_XS);
	data_CFD_Ys.resize(species_set.NS, YsNameMap); data_CFD_Ys.asgName(VAR_VECTOR_YS);



	/*
	 * Assign sample data format for
	 * 	- diffusion coefficient
	 * 	- viscosity coefficient
	 */
	Common::Map1D<std::string, int>	munameMap (species_set.NS);
	Common::Map1D<std::string, int>	DnameMap (species_set.NS);

	for (int s = 0; s <= species_set.NS-1; s++)
	{
		munameMap.insert("mu_" + species_set.species[s].name, s);
		DnameMap.insert("D_" + species_set.species[s].name, s);
	}

	data_CFD_diffusion_coeff.resize(species_set.NS, DnameMap); data_CFD_diffusion_coeff.asgName("Species diffusion coefficient");
	data_CFD_viscosity_coeff.resize(species_set.NS, DnameMap); data_CFD_viscosity_coeff.asgName("Species viscosity coefficient");


	/*
	 * Assign Sample data format for
	 *	 - dp/dQ
	 */
	Common::Map1D<std::string, int>	dpdQMap (numData);
	for (int i = 0; i <= numData-1; i++)	dpdQMap.insert("dp/"+Q[i], i);
	data_CFD_dp_dQ.resize(numData, dpdQMap); data_CFD_dp_dQ.asgName("dp/dQ");


	/*
	 * Assign sample data format for
	 * 	- div Vc
	 */
	Common::Map1D<std::string, int>	divVcMap (1);
	divVcMap.insert("divergence V", 0);
	data_CFD_divVc.resize(1, divVcMap); data_CFD_divVc.asgName("divergence Vc");


	// 2D Data
	/*
	 * Assign Sample data format for
	 * 	 - Jacobian_inviscid_plus
	 *	 - Jacobian_inviscid_minus
	 *
	 *	 - Jacobian_viscous_plus
	 *	 - Jacobian_viscous_minus
	 *
	 *	 - Jacobian_source
	 *
	 *	 - dT/dQ
	 */
	vector<string> T(ne);
	vector<string> kappa(ne);

	ne	= 0;
	T[ne]		= "T";
	kappa[ne]	= "k";
	ne++;

	if (problem_setup.NER != 0)
	{
		T[ne]		= "T_rot";
		kappa[ne]	= "k_rot";
		ne++;
	}

	if (problem_setup.NEV != 0)
	{
		T[ne]		= "T_vib";
		kappa[ne]	= "k_vib";
		ne++;
	}

	if (problem_setup.NEE != 0)
	{
		T[ne]		= "T_e";
		kappa[ne]	= "k_e";
		ne++;
	}

	Common::Map2D<std::string, string, int>	JacobianNameMap (numData*numData);
	Common::Map2D<std::string, string, int>	JacobianSourceNameMap (numData*numData);
	Common::Map2D<std::string, string, int>	dTdQNameMap (ne*numData);
	Common::Map2D<std::string, string, int>	kappaNameMap (ne*species_set.NS);

	int ii = 0;
	for (int j = 0; j <= numData-1; j++)
	{
		for (int i = 0; i <= numData-1; i++)
		{
			JacobianNameMap.insert("d"+F[i], "d"+Q[j], ii);
			JacobianSourceNameMap.insert("dS_"+Q[i], "d"+Q[j], ii);
			ii++;
		}
	}

	data_CFD_Jacobian_inviscid_minus.resize(numData, numData, JacobianNameMap); data_CFD_Jacobian_inviscid_minus.asgName("Inviscid Jacobian minus");
	data_CFD_Jacobian_inviscid_plus.resize(numData, numData, JacobianNameMap); data_CFD_Jacobian_inviscid_plus.asgName("Inviscid Jacobian plus");

	data_CFD_Jacobian_viscous_minus.resize(numData, numData, JacobianNameMap); data_CFD_Jacobian_viscous_minus.asgName("Viscous Jacobian minus");
	data_CFD_Jacobian_viscous_plus.resize(numData, numData, JacobianNameMap); data_CFD_Jacobian_viscous_plus.asgName("Viscous Jacobian plus");

	data_CFD_Jacobian_source.resize(numData, numData, JacobianSourceNameMap); data_CFD_Jacobian_source.asgName("Source term Jacobian");

	ii = 0;
	for (int j = 0; j <= numData-1; j++)
	{
		for (int i = 0; i <= ne-1; i++)
		{
			dTdQNameMap.insert("d"+T[i], "d"+Q[j], ii);
			ii++;
		}
	}
	data_CFD_dT_dQ.resize(ne, numData, dTdQNameMap); data_CFD_dT_dQ.asgName("dT / dQ");

	ii = 0;
	for (int j = 0; j <= species_set.NS-1; j++)
	{
		for (int i = 0; i <= ne-1; i++)
		{
			kappaNameMap.insert(kappa[i], species_set.species[j].name, ii);
			ii++;
		}
	}
	data_CFD_thermal_conductivity_coeff.resize(ne, species_set.NS, kappaNameMap); data_CFD_thermal_conductivity_coeff.asgName("Species thermal conductivity");
}
