/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 15, 2015
 *      			Author: Minkwan Kim
 *
 * AppliationConstants.hpp
 * 			-  
 *  
 */
#ifndef APPLIATIONCONSTANTS_HPP_
#define APPLIATIONCONSTANTS_HPP_


#define VAR_VECTOR_Q	"Conservative Variables"
#define VAR_VECTOR_V	"Primitive Variables"
#define VAR_VECTOR_W	"W"
#define VAR_VECTOR_MIX	"Mixture variables"
#define VAR_VECTOR_XS	"Species mole fraction"
#define VAR_VECTOR_YS	"Species mass fraction"
#define VAR_VECTOR_F_INV	"Inviscid Flux"
#define VAR_VECTOR_F_VIS	"Viscous Flux"


/*
 * Q/V/W variables names
 */
#define VAR_RHO_PRE		"rho<sub>"
#define VAR_RHO_POST	"</sub> [kg/m<sup>3</sup>]"

#define VAR_RHO_U		"rho_u"
#define VAR_RHO_V		"rho_v"
#define VAR_RHO_W		"rho_w"

#define VAR_U			"u [m/s]"
#define VAR_V			"v [m/s]"
#define VAR_W			"w [m/s]"

#define VAR_T			"T <sub>tra</sub> [K]"
#define VAR_T_ROT		"T <sub>rot</sub> [K]"
#define VAR_T_VIB		"T <sub>vib</sub> [K]"
#define VAR_T_E			"T <sub>e</sub> [K]"
#define VAR_P			"P [N/m<sup>2</sup>]"

#define VAR_E			"E <sub>total</sub> [J]"
#define VAR_E_ROT		"E <sub>rot</sub> [J]"
#define VAR_E_VIB		"E <sub>vib</sub> [J]"
#define VAR_E_E			"E <sub>e</sub> [J]"


#define VAR_FLUX		"F_"
#define VAR_F_INV		"Finv_"
#define VAR_F_VIS		"Fvis_"

#define VAR_RHO			"rho_"
#define VAR_SPECIES		"s"
#define VAR_U_WO_UNIT	"u"
#define VAR_V_WO_UNIT	"v"
#define VAR_W_WO_UNIT	"w"

#define VAR_T_WO_UNIT		"T"
#define VAR_T_ROT_WO_UNIT	"T_rot"
#define VAR_T_VIB_WO_UNIT	"T_vib"
#define VAR_T_E_WO_UNIT		"T_e"

#define VAR_E_WO_UNIT		"E"
#define VAR_E_ROT_WO_UNIT	"E_rot"
#define VAR_E_VIB_WO_UNIT	"E_vib"
#define VAR_E_E_WO_UNIT		"E_e"

/*
 * Mixture variables
 */
#define	VAR_MIX_SUB			"_mix"
#define	VAR_TRA_SUB			"_tra"
#define	VAR_ROT_SUB			"_rot"
#define	VAR_VIB_SUB			"_vib"
#define	VAR_ELE_SUB			"_ele"
#define	VAR_E_SUB			"_e"


#define	VAR_MIX				"mix"
#define VAR_R				"R"
#define VAR_M				"M"
#define VAR_GAMMA			"gamma"
#define VAR_CV_TRA			"CV_tra"
#define VAR_CV_ROT			"CV_rot"
#define VAR_CV_VIB			"CV_vib"
#define VAR_CV_ELE			"CV_ele"
#define VAR_CV_TR			"CV_tr"

#define VAR_VISCOSITY_COEFF	"mu"
#define VAR_THERMAL_COND	"kappa"
#define VAR_VISCOSITY_COEFF	"mu"
#define VAR_DIFFUSION_COEFF	"D"


#define VAR_XS		"Xs_"
#define VAR_YS		"Ys_"







#endif /* APPLIATIONCONSTANTS_HPP_ */
