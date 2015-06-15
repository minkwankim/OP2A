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


	for (int s = 0; s <= species_set.NS-1; s++)
	{
		Q[s]	= "rho_" + species_set.species[s].name;
		V[s]	= "rho_" + species_set.species[s].name;
		F[s]	= "F_rho_" + species_set.species[s].name;

		QnameMap.insert("rho<sub>" + species_set.species[s].name + "</sub> [kg/m<sup>3</sup>]", s);
		VnameMap.insert("rho<sub>" + species_set.species[s].name + "</sub> [kg/m<sup>3</sup>]", s);
		WnameMap.insert("rho<sub>" + species_set.species[s].name + "</sub> [kg/m<sup>3</sup>]", s);

		FluxInvnameMap.insert("Finv_s" + Common::StringOps::to_str<int>(s), s);
		FluxVisnameMap.insert("Fvis_s" + Common::StringOps::to_str<int>(s), s);
	}

	Q[species_set.NS]	= "rho_u";
	V[species_set.NS]	= "u";
	F[species_set.NS]	= "F_rho_u";
	QnameMap.insert("rho u", species_set.NS);
	VnameMap.insert("u [m/s]", species_set.NS);
	WnameMap.insert("u [m/s]", species_set.NS);

	FluxInvnameMap.insert("Finv_rhou", species_set.NS);
	FluxVisnameMap.insert("Fvis_rhou", species_set.NS);


	Q[species_set.NS+1]	= "rho_v";
	V[species_set.NS+1]	= "v";
	F[species_set.NS+1]	= "F_rho_v";
	QnameMap.insert("rho v", species_set.NS+1);
	VnameMap.insert("v [m/s]", species_set.NS+1);
	WnameMap.insert("v [m/s]", species_set.NS+1);

	FluxInvnameMap.insert("Finv_rhov", species_set.NS+1);
	FluxVisnameMap.insert("Fvis_rhov", species_set.NS+1);


	if (grid.ND == 3)
	{
		Q[species_set.NS+2]	= "rho_w";
		V[species_set.NS+2]	= "w";
		F[species_set.NS+2]	= "F_rho_w";

		QnameMap.insert("rho w", species_set.NS+2);
		VnameMap.insert("w [m/s]", species_set.NS+2);
		WnameMap.insert("w [m/s]", species_set.NS+2);

		FluxInvnameMap.insert("Finv_rhow", species_set.NS+2);
		FluxVisnameMap.insert("Fvis_rhow", species_set.NS+2);
	}

	Q[species_set.NS+grid.ND]	= "E";
	V[species_set.NS+grid.ND]	= "T";
	F[species_set.NS+grid.ND]	= "F_E";
	QnameMap.insert("E <sub>total</sub> [J]", species_set.NS+grid.ND);
	VnameMap.insert("T <sub>tra</sub> [K]", species_set.NS+grid.ND);
	WnameMap.insert("P [N/m<sup>2</sup>]", species_set.NS+grid.ND);


	int ne	= 0;
	FluxInvnameMap.insert("Finv_E", species_set.NS+grid.ND+ne);
	FluxVisnameMap.insert("Fvis_E", species_set.NS+grid.ND+ne);
	ne++;

	if (problem_setup.NER != 0)
	{
		Q[species_set.NS+grid.ND+ne]	= "E_rot";
		V[species_set.NS+grid.ND+ne]	= "T_rot";
		F[species_set.NS+grid.ND+ne]	= "F_Erot";
		QnameMap.insert("E <sub>rot</sub> [J]", species_set.NS+grid.ND+ne);
		VnameMap.insert("T <sub>rot</sub> [K]", species_set.NS+grid.ND+ne);
		WnameMap.insert("T <sub>rot</sub> [K]", species_set.NS+grid.ND+ne);

		FluxInvnameMap.insert("Finv_Erot", species_set.NS+grid.ND+ne);
		FluxVisnameMap.insert("Fvis_Erot", species_set.NS+grid.ND+ne);

		ne++;
	}

	if (problem_setup.NEV != 0)
	{
		Q[species_set.NS+grid.ND+ne]	= "E_vib";
		V[species_set.NS+grid.ND+ne]	= "T_vib";
		F[species_set.NS+grid.ND+ne]	= "F_Evib";
		QnameMap.insert("E <sub>vib</sub> [J]", species_set.NS+grid.ND+ne);
		VnameMap.insert("T <sub>vib</sub> [K]", species_set.NS+grid.ND+ne);
		WnameMap.insert("T <sub>vib</sub> [K]", species_set.NS+grid.ND+ne);

		FluxInvnameMap.insert("Finv_Evib", species_set.NS+grid.ND+ne);
		FluxVisnameMap.insert("Fvis_Evib", species_set.NS+grid.ND+ne);

		ne++;
	}

	if (problem_setup.NEE != 0)
	{
		Q[species_set.NS+grid.ND+ne]	= "E_e";
		V[species_set.NS+grid.ND+ne]	= "T_e";
		F[species_set.NS+grid.ND+ne]	= "F_Ee";
		QnameMap.insert("E <sub>e</sub> [J]", species_set.NS+grid.ND+ne);
		VnameMap.insert("T <sub>e</sub> [K]", species_set.NS+grid.ND+ne);
		WnameMap.insert("T <sub>e</sub> [K]", species_set.NS+grid.ND+ne);

		FluxInvnameMap.insert("Finv_Ee", species_set.NS+grid.ND+ne);
		FluxVisnameMap.insert("Fvis_Ee", species_set.NS+grid.ND+ne);

		ne++;
	}


	data_CFD_Q.resize(numData, QnameMap);	data_CFD_Q.asgName("Conservative variables");
	data_CFD_V.resize(numData, VnameMap);	data_CFD_V.asgName("Primitive variables");
	data_CFD_W.resize(numData, WnameMap);	data_CFD_W.asgName("W");

	data_CFD_Flux_inviscid.resize(numData, FluxInvnameMap);	data_CFD_Flux_inviscid.asgName("Inviscid Flux");
	data_CFD_Flux_viscous.resize(numData, FluxVisnameMap);	data_CFD_Flux_viscous.asgName("Viscous Flux");



	/*
	 * Assign Sample data format for
	 * 	 - Mixture data
	 */

	int nCv	= 1;
	if (problem_setup.NER != 0)	nCv++;

	Common::Map1D<std::string, int>	MixtureNameMap (4+nCv+1+ne);

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

	MixtureNameMap.insert("mu_mix", 4+nCv);

	ne	= 0;

	MixtureNameMap.insert("kappa_tra", 4+nCv+ne);
	ne++;

	if (problem_setup.NER != 0)
	{
		MixtureNameMap.insert("kappa_rot", 4+nCv+ne);
		ne++;
	}

	if (problem_setup.NEV != 0)
	{
		MixtureNameMap.insert("kappa_vib", 4+nCv+ne);
		ne++;
	}

	if (problem_setup.NEE != 0)
	{
		MixtureNameMap.insert("kappa_ele", 4+nCv+ne);
		ne++;
	}

	data_CFD_mixture.resize(4+nCv+1+ne, MixtureNameMap); data_CFD_Q.asgName("Mixture variables");



	/*
	 * Assign Sample data format for
	 * 	 - Xs
	 * 	 - Ys
	 */
	Common::Map1D<std::string, int>	YsNameMap (species_set.NS);
	Common::Map1D<std::string, int>	XsNameMap (species_set.NS);

	for (int s = 0; s <= species_set.NS-1; s++)
	{
		YsNameMap.insert("Y<sub>" + species_set.species[s].name + "</sub>", s);
		XsNameMap.insert("X<sub>" + species_set.species[s].name + "</sub>", s);
	}

	data_CFD_Xs.resize(species_set.NS, XsNameMap); data_CFD_Xs.asgName("Species mole fraction");
	data_CFD_Ys.resize(species_set.NS, YsNameMap); data_CFD_Ys.asgName("Species mass fraction");



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
