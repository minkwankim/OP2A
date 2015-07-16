/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jul 7, 2015
 *      			Author: Minkwan Kim
 *
 * DerivativesCaseA.cpp
 * 			-  
 *  
 */


#include <omp.h>

#include "Common/include/MultiDimension.hpp"
#include "CFD/include/Derivatives.hpp"


namespace OP2A{
namespace CFD{


void Derivatives::d2pdQ_CaseA1(Data::DataStorageVector<Data::DataStorage>& data1D,
								Data::DataStorage& dT,
								CHEM::SpeciesSet& species_set, int ND,
								unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX,
								Data::DataStorage2D& d2p)
{
	int VAR	= data1D(indexQ).numData;
	double rho		= data1D(indexMIX)(0);
	double Cv_bar	= data1D(indexMIX)(4) / rho;
	double R_bar	= data1D(indexMIX)(1) / rho;
	double R_bar_over_Cv_bar	= R_bar / Cv_bar;

	double u_square	= 0.0;
	for (int k = 0; k <= ND-1; k++)	u_square	+= pow(data1D(indexV)(species_set.NS+k), 2.0);

	// Variable 1
	vector<double>	dR_o_Cv_dQ(VAR, 0.0);
	double rho_Cv_bar2	= rho * pow(Cv_bar, 2.0);

#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		dR_o_Cv_dQ[s]	= (species_set.species[s].R*Cv_bar - species_set.species[s].Cv_tr*R_bar) / rho_Cv_bar2;
	}


	// Variable 2
	vector<double>	u_square_dQ(VAR, 0.0);
	double temp1	= -2.0*u_square/rho;

#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		u_square_dQ[s] = temp1;
	}

#pragma ivdep
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		u_square_dQ[k]	= 2.0*data1D(indexV)(k)/rho;
	}

	vector <vector<double> >	du_dQ	= Common::vector_2D(ND, VAR, 0.0);
	for (int k = 0; k <= ND-1; k++)
	{
		temp1	= -data1D(indexV)(species_set.NS+k)/rho;

#pragma ivdep
		for (int s = 0; s <= species_set.NS-1; s++)
		{
			du_dQ[k][s] = temp1;
		}

		du_dQ[k][species_set.NS+k]	= 1.0/rho;
	}



	/*
	 * Calculate d2p_dQ2
	 */
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		temp1	= 0.5*u_square - species_set.species[s].h0;

		double temp2 = species_set.species[s].R - R_bar_over_Cv_bar*species_set.species[s].Cv_tr;
		double temp3 = species_set.species[s].Cv_tr*data1D(indexV)(species_set.NS+ND);

#pragma ivdep
		for (int k = 0; k <= VAR-1; k++)
		{
			d2p(s,k) = dR_o_Cv_dQ[k]*temp1 + 0.5*R_bar_over_Cv_bar*u_square_dQ[k];
			d2p(s,k) += temp2 * dT(k) -  dR_o_Cv_dQ[k]*temp3;
		}
	}

	for (int k = 0; k <= ND-1; k++)
	{

#pragma ivdep
		for (int i = 0; i <= VAR-1; i++)
		{
			d2p(species_set.NS+k,i)	= -(dR_o_Cv_dQ[i]*data1D(indexV)(species_set.NS+k) + R_bar_over_Cv_bar*du_dQ[k][i]);
		}
	}

#pragma ivdep
	for (int i = 0; i <= VAR-1; i++)
	{
		d2p(species_set.NS+ND, i)	= dR_o_Cv_dQ[i];
	}

	for (int k = species_set.NS+ND+1; k <= VAR-1; k++)
	{
#pragma ivdep
		for (int i = 0; i <= VAR-1; i++)
		{
			d2p(k,i) = -dR_o_Cv_dQ[i];
		}
	}
}


void Derivatives::d2pdQ_CaseA2(Data::DataStorageVector<Data::DataStorage>& data1D,
								Data::DataStorage& dT,
								CHEM::SpeciesSet& species_set, int ND,
								unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX,
								Data::DataStorage2D& d2p)
{
	int VAR	= data1D(indexQ).numData;

	double rho		= data1D(indexMIX)(0);
	double Cv_bar	= data1D(indexMIX)(4) / rho;
	double R_bar	= data1D(indexMIX)(1) / rho;
	double R_bar_over_Cv_bar	= R_bar / Cv_bar;

	double u_square	= 0.0;
	for (int k = 0; k <= ND-1; k++)	u_square	+= pow(data1D(indexV)(species_set.NS+k), 2.0);

	// Variable 1
	vector<double>	dR_o_Cv_dQ(VAR, 0.0);
	double rho_Cv_bar2	= rho * pow(Cv_bar, 2.0);

#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		dR_o_Cv_dQ[s]	= (species_set.species[s].R*Cv_bar - species_set.species[s].Cv_tra*R_bar) / rho_Cv_bar2;
	}


	// Variable 2
	vector<double>	u_square_dQ(VAR, 0.0);
	double temp1	= -2.0*u_square/rho;

#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		u_square_dQ[s] = temp1;
	}

#pragma ivdep
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		u_square_dQ[k]	= 2.0*data1D(indexV)(k)/rho;
	}

	vector <vector<double> >	du_dQ	= Common::vector_2D(ND, VAR, 0.0);
	for (int k = 0; k <= ND-1; k++)
	{
		temp1	= -data1D(indexV)(species_set.NS+k)/rho;

#pragma ivdep
		for (int s = 0; s <= species_set.NS-1; s++)
		{
			du_dQ[k][s] = temp1;
		}

		du_dQ[k][species_set.NS+k]	= 1.0/rho;
	}



	/*
	 * Calculate d2p_dQ2
	 */
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		temp1	= 0.5*u_square - species_set.species[s].h0;

		double temp2 = species_set.species[s].R - R_bar_over_Cv_bar*species_set.species[s].Cv_tra;
		double temp3 = species_set.species[s].Cv_tra*data1D(indexV)(species_set.NS+ND);

#pragma ivdep
		for (int k = 0; k <= VAR-1; k++)
		{
			d2p(s,k) = dR_o_Cv_dQ[k]*temp1 + 0.5*R_bar_over_Cv_bar*u_square_dQ[k];
			d2p(s,k) += temp2 * dT(k) -  dR_o_Cv_dQ[k]*temp3;
		}
	}


	for (int k = 0; k <= ND-1; k++)
	{

#pragma ivdep
		for (int i = 0; i <= VAR-1; i++)
		{
			d2p(species_set.NS+k,i)	= -(dR_o_Cv_dQ[i]*data1D(indexV)(species_set.NS+k) + R_bar_over_Cv_bar*du_dQ[k][i]);
		}
	}


#pragma ivdep
	for (int i = 0; i <= VAR-1; i++)
	{
		d2p(species_set.NS+ND, i)	= dR_o_Cv_dQ[i];
	}


	for (int k = species_set.NS+ND+1; k <= VAR-1; k++)
	{
#pragma ivdep
		for (int i = 0; i <= VAR-1; i++)
		{
			d2p(k,i) = -dR_o_Cv_dQ[i];
		}
	}
}







void Derivatives::d2pdQ_CaseA(Data::DataStorageVector<Data::DataStorage>& data1D,
									Data::DataStorage& dT,
									Data::DataStorage2D& d2T,
									CHEM::SpeciesSet& species_set, int ND,
									unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX,
									Data::DataStorage2D& d2p)
{
	int VAR	= data1D(indexQ).numData;
	double rhoR	= data1D(indexMIX)(1);
	vector<double>	d_rhoR (VAR, 0.0);

#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)	d_rhoR[s]	= species_set.species[s].R;

#pragma ivdep
	for (int j = 0; j <= VAR-1; j++)
	{
		for (int i = 0; i <= VAR-1; i++)
		{
			d2p(i,j) = d_rhoR[i]*dT(j) + rhoR*d2T(i,j) + d_rhoR[j]*dT(i);
		}
	}
}






void Derivatives::d2pdQ_CaseB1(Data::DataStorageVector<Data::DataStorage>& data1D,
								Data::DataStorage& dT,
								Data::DataStorage& dTe,
								CHEM::SpeciesSet& species_set, int ND,
								unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX,
								Data::DataStorage2D& d2p)
{
	int VAR	= data1D(indexQ).numData;

	double rho			= data1D(indexMIX)(0);
	double Cv_bar_wo_e	= data1D(indexMIX)(4);
	double R_bar_wo_e	= data1D(indexMIX)(1);
	double R_bar_over_Cv_bar	= R_bar_wo_e / Cv_bar_wo_e;

	double u_square	= 0.0;
	for (int k = 0; k <= ND-1; k++)	u_square	+= pow(data1D(indexV)(species_set.NS+k), 2.0);


	// Variable 1
	vector<double>	dR_o_Cv_dQ(VAR, 0.0);
	double rho_Cv_bar2	= pow(Cv_bar_wo_e, 2.0);
	double rhoe_Re;
	double rhoe_Cve;

#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			dR_o_Cv_dQ[s]	= (species_set.species[s].R*Cv_bar_wo_e - species_set.species[s].Cv_tr*R_bar_wo_e) / rho_Cv_bar2;
		}
		else
		{
			dR_o_Cv_dQ[s]	= 0.0;
			rhoe_Re			= data1D(indexQ)(s)*species_set.species[s].R;
			rhoe_Cve		= data1D(indexQ)(s)*species_set.species[s].Cv_tra;
		}
	}


	// Variable 2
	vector<double>	u_square_dQ(VAR, 0.0);
	double temp1	= -2.0*u_square/rho;

#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		u_square_dQ[s] = temp1;
	}

#pragma ivdep
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		u_square_dQ[k]	= 2.0*data1D(indexV)(k)/rho;
	}

	vector <vector<double> >	du_dQ	= Common::vector_2D(ND, VAR, 0.0);
	for (int k = 0; k <= ND-1; k++)
	{
		temp1	= -data1D(indexV)(species_set.NS+k)/rho;

#pragma ivdep
		for (int s = 0; s <= species_set.NS-1; s++)
		{
			du_dQ[k][s] = temp1;
		}

		du_dQ[k][species_set.NS+k]	= 1.0/rho;
	}



	/*
	 * Calculate d2p_dQ2
	 */
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			temp1	= 0.5*u_square - species_set.species[s].h0;

			double temp2 = species_set.species[s].R - R_bar_over_Cv_bar*species_set.species[s].Cv_tr;
			double temp3 = species_set.species[s].Cv_tr*data1D(indexV)(species_set.NS+ND);

	#pragma ivdep
			for (int k = 0; k <= VAR-1; k++)
			{
				d2p(s,k) = dR_o_Cv_dQ[k]*temp1 + 0.5*R_bar_over_Cv_bar*u_square_dQ[k];
				d2p(s,k) += temp2 * dT(k) -  dR_o_Cv_dQ[k]*temp3;
			}
		}
		else
		{
			temp1	= 0.5*u_square;

	#pragma ivdep
			for (int k = 0; k <= VAR-1; k++)
			{
				d2p(s,k) = dR_o_Cv_dQ[k]*temp1 + 0.5*R_bar_over_Cv_bar*u_square_dQ[k];
			}
		}
	}

	for (int k = 0; k <= ND-1; k++)
	{

#pragma ivdep
		for (int i = 0; i <= VAR-1; i++)
		{
			d2p(species_set.NS+k,i)	= -(dR_o_Cv_dQ[i]*data1D(indexV)(species_set.NS+k) + R_bar_over_Cv_bar*du_dQ[k][i]);
		}
	}


	// Total Energy
#pragma ivdep
	for (int i = 0; i <= VAR-1; i++)
	{
		d2p(species_set.NS+ND, i)	= dR_o_Cv_dQ[i];
	}


	// Other energy mode
	for (int k = species_set.NS+ND+1; k <= VAR-1; k++)
	{
#pragma ivdep
		for (int i = 0; i <= VAR-1; i++)
		{
			d2p(k,i) = -dR_o_Cv_dQ[i];
		}
	}
}


void Derivatives::d2pdQ_CaseB2(Data::DataStorageVector<Data::DataStorage>& data1D,
								Data::DataStorage& dT,
								Data::DataStorage& dTe,
								CHEM::SpeciesSet& species_set, int ND,
								unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX,
								Data::DataStorage2D& d2p)
{
	int VAR	= data1D(indexQ).numData;

	double rho			= data1D(indexMIX)(0);
	double Cv_bar_wo_e	= data1D(indexMIX)(4);
	double R_bar_wo_e	= data1D(indexMIX)(1);
	double R_bar_over_Cv_bar	= R_bar_wo_e / Cv_bar_wo_e;

	double u_square	= 0.0;
	for (int k = 0; k <= ND-1; k++)	u_square	+= pow(data1D(indexV)(species_set.NS+k), 2.0);


	// Variable 1
	vector<double>	dR_o_Cv_dQ(VAR, 0.0);
	double rho_Cv_bar2	= pow(Cv_bar_wo_e, 2.0);
	double rhoe_Re;
	double rhoe_Cve;

#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			dR_o_Cv_dQ[s]	= (species_set.species[s].R*Cv_bar_wo_e - species_set.species[s].Cv_tra*R_bar_wo_e) / rho_Cv_bar2;
		}
		else
		{
			dR_o_Cv_dQ[s]	= 0.0;
			rhoe_Re			= data1D(indexQ)(s)*species_set.species[s].R;
			rhoe_Cve		= data1D(indexQ)(s)*species_set.species[s].Cv_tra;
		}
	}


	// Variable 2
	vector<double>	u_square_dQ(VAR, 0.0);
	double temp1	= -2.0*u_square/rho;

#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		u_square_dQ[s] = temp1;
	}

#pragma ivdep
	for (int k = species_set.NS; k <= species_set.NS+ND-1; k++)
	{
		u_square_dQ[k]	= 2.0*data1D(indexV)(k)/rho;
	}

	vector <vector<double> >	du_dQ	= Common::vector_2D(ND, VAR, 0.0);
	for (int k = 0; k <= ND-1; k++)
	{
		temp1	= -data1D(indexV)(species_set.NS+k)/rho;

#pragma ivdep
		for (int s = 0; s <= species_set.NS-1; s++)
		{
			du_dQ[k][s] = temp1;
		}

		du_dQ[k][species_set.NS+k]	= 1.0/rho;
	}



	/*
	 * Calculate d2p_dQ2
	 */
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if(species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			temp1	= 0.5*u_square - species_set.species[s].h0;

			double temp2 = species_set.species[s].R - R_bar_over_Cv_bar*species_set.species[s].Cv_tra;
			double temp3 = species_set.species[s].Cv_tra*data1D(indexV)(species_set.NS+ND);

	#pragma ivdep
			for (int k = 0; k <= VAR-1; k++)
			{
				d2p(s,k) = dR_o_Cv_dQ[k]*temp1 + 0.5*R_bar_over_Cv_bar*u_square_dQ[k];
				d2p(s,k) += temp2 * dT(k) -  dR_o_Cv_dQ[k]*temp3;
			}
		}
		else
		{
			temp1	= 0.5*u_square;

	#pragma ivdep
			for (int k = 0; k <= VAR-1; k++)
			{
				d2p(s,k) = dR_o_Cv_dQ[k]*temp1 + 0.5*R_bar_over_Cv_bar*u_square_dQ[k];
			}
		}
	}

	for (int k = 0; k <= ND-1; k++)
	{

#pragma ivdep
		for (int i = 0; i <= VAR-1; i++)
		{
			d2p(species_set.NS+k,i)	= -(dR_o_Cv_dQ[i]*data1D(indexV)(species_set.NS+k) + R_bar_over_Cv_bar*du_dQ[k][i]);
		}
	}


	// Total Energy
#pragma ivdep
	for (int i = 0; i <= VAR-1; i++)
	{
		d2p(species_set.NS+ND, i)	= dR_o_Cv_dQ[i];
	}


	// Other energy mode
	for (int k = species_set.NS+ND+1; k <= VAR-1; k++)
	{
#pragma ivdep
		for (int i = 0; i <= VAR-1; i++)
		{
			d2p(k,i) = -dR_o_Cv_dQ[i];
		}
	}
}








void Derivatives::d2pdQ_CaseB(Data::DataStorageVector<Data::DataStorage>& data1D,
								Data::DataStorage& dT,
								Data::DataStorage& dTe,
								Data::DataStorage2D& d2T,
								Data::DataStorage2D& d2Te,
								CHEM::SpeciesSet& species_set, int ND,
								unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX,
								Data::DataStorage2D& d2p)
{
	int VAR	= data1D(indexQ).numData;

	double rhoR	= data1D(indexMIX)(1);
	double rhoeRe;

	vector<double>	d_rhoR (VAR, 0.0);
	vector<double>	d_rhoeRe (VAR, 0.0);

#pragma ivdep
	for (int s = 0; s <= species_set.NS-1; s++)
	{
		if (species_set.species[s].type != CHEM::SpeciesType::Electron)
		{
			d_rhoR[s]	= species_set.species[s].R;
		}
		else
		{
			d_rhoeRe[s]	= species_set.species[s].R;
			rhoeRe	= data1D(indexQ)(s) * species_set.species[s].R;
		}
	}

#pragma ivdep
	for (int j = 0; j <= VAR-1; j++)
	{
		for (int i = 0; i <= VAR-1; i++)
		{
			d2p(i,j) = d_rhoR[i]*dT(j) + rhoR*d2T(i,j) + d_rhoR[j]*dT(i)
					+ d_rhoeRe[i]*dTe(j) + rhoeRe*d2Te(i,j) + d_rhoeRe[j]*dTe(j);
		}
	}
}



void Derivatives::d2pdQ(Data::DataStorageVector<Data::DataStorage>& data1D,
				  Data::DataStorageVector<Data::DataStorage>& dT,
				  CHEM::SpeciesSet& species_set, int ND,
				  unsigned int type, unsigned int indexQ, unsigned int indexV, unsigned int indexW,  unsigned int indexMIX,
				  Data::DataStorage2D& d2p)
{
	int indexT  = 0;
	int indexTe = 0;
	int VAR = data1D(indexQ).numData;

	Data::DataStorage2D d2T(VAR, VAR);
	Data::DataStorage2D d2Te(VAR, VAR);

	switch (type)
	{
	case 1:
		//DerivativesType1::d2TdQ2(data1D(indexQ), data1D(indexV), data1D(indexMIX), species_set, ND, dT(indexT), d2T);
		//d2pdQ_CaseA(data1D, dT(indexT), d2T, species_set, ND, indexQ, indexV, indexW,  indexMIX, d2p);
		d2pdQ_CaseA1(data1D, dT(indexT), species_set, ND, indexQ, indexV, indexW,  indexMIX, d2p);
		break;

	case 2:
		indexTe = 1;
		//DerivativesType2::d2TdQ2(data1D(indexQ), data1D(indexV), data1D(indexMIX), species_set, ND, dT(indexT), dT(indexTe), d2T, d2Te);
		//d2pdQ_CaseB(data1D, dT(indexT), dT(indexTe), d2T, d2Te, species_set, ND, indexQ, indexV, indexW,  indexMIX, d2p);
		d2pdQ_CaseB1(data1D, dT(indexT), dT(indexTe), species_set, ND, indexQ, indexV, indexW,  indexMIX, d2p);
		break;

	case 3:
		//DerivativesType3::d2TdQ2(data1D(indexQ), data1D(indexV), data1D(indexMIX), species_set, ND, dT(indexT), d2T);
		//d2pdQ_CaseA(data1D, dT(indexT), d2T, species_set, ND, indexQ, indexV, indexW,  indexMIX, d2p);
		d2pdQ_CaseA1(data1D, dT(indexT), species_set, ND, indexQ, indexV, indexW,  indexMIX, d2p);
		break;

	case 4:
		indexTe = 2;
		//DerivativesType4::d2TdQ2(data1D(indexQ), data1D(indexV), data1D(indexMIX), species_set, ND, dT(indexT), dT(indexTe), d2T, d2Te);
		//d2pdQ_CaseB(data1D, dT(indexT), dT(indexTe), d2T, d2Te, species_set, ND, indexQ, indexV, indexW,  indexMIX, d2p);
		d2pdQ_CaseB1(data1D, dT(indexT), dT(indexTe), species_set, ND, indexQ, indexV, indexW,  indexMIX, d2p);
		break;

	case 5:
		//DerivativesType5::d2TdQ2(data1D(indexQ), data1D(indexV), data1D(indexMIX), species_set, ND, dT(indexT), d2T);
		//d2pdQ_CaseA(data1D, dT(indexT), d2T, species_set, ND, indexQ, indexV, indexW,  indexMIX, d2p);
		d2pdQ_CaseA2(data1D, dT(indexT), species_set, ND, indexQ, indexV, indexW,  indexMIX, d2p);
		break;

	case 6:
		indexTe = 2;
		//DerivativesType6::d2TdQ2(data1D(indexQ), data1D(indexV), data1D(indexMIX), species_set, ND, dT(indexT), dT(indexTe), d2T, d2Te);
		//d2pdQ_CaseB(data1D, dT(indexT), dT(indexTe), d2T, d2Te, species_set, ND, indexQ, indexV, indexW,  indexMIX, d2p);
		d2pdQ_CaseB2(data1D, dT(indexT), dT(indexTe), species_set, ND, indexQ, indexV, indexW,  indexMIX, d2p);
		break;

	case 7:
		//DerivativesType7::d2TdQ2(data1D(indexQ), data1D(indexV), data1D(indexMIX), species_set, ND, dT(indexT), d2T);
		//d2pdQ_CaseA(data1D, dT(indexT), d2T, species_set, ND, indexQ, indexV, indexW,  indexMIX, d2p);
		d2pdQ_CaseA2(data1D, dT(indexT), species_set, ND, indexQ, indexV, indexW,  indexMIX, d2p);
		break;

	case 8:
		indexTe = 3;
		//DerivativesType8::d2TdQ2(data1D(indexQ), data1D(indexV), data1D(indexMIX), species_set, ND, dT(indexT), dT(indexTe), d2T, d2Te);
		//d2pdQ_CaseB(data1D, dT(indexT), dT(indexTe), d2T, d2Te, species_set, ND, indexQ, indexV, indexW,  indexMIX, d2p);
		d2pdQ_CaseB2(data1D, dT(indexT), dT(indexTe), species_set, ND, indexQ, indexV, indexW,  indexMIX, d2p);
		break;
	}
}




}
}
