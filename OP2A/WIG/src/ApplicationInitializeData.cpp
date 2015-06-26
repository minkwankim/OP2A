/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 15, 2015
 *      			Author: Minkwan Kim
 *
 * ApplicationInitializeData.cpp
 * 			-  
 *  
 */


#include "Common/include/CodeLocation.hpp"
#include "Common/include/Exception_NPExceed.hpp"
#include "Common/include/Exception_DimensionMatch.hpp"

#include "Common/include/StringOps.hpp"
#include "Common/include/DataReadFromFile.hpp"

#include "CFD/include/VariableSamplesWIG.hpp"

#include "CFD/include/VariableConstants.hpp"
#include "CFD/include/VariableChange.hpp"

#include "../include/OP2A_Application.hpp"




void ApplicationOP2A::InitializeData(unsigned int num_ic, bool use_restart_file)
{
	CalculateIC();


	int indexQ = grid.cells[1].data1D.dataMap.find(NAME_Q);
	int indexV = grid.cells[1].data1D.dataMap.find(NAME_V);

	if (use_restart_file == true)
	{
		cout << " ---> Initializing using restart file......" << endl;
		string line;
		string restart_filename;

		restart_filename	= problem_setup.name + ".rst";
		ifstream restart_file(restart_filename.c_str());

		if (!restart_file)
		{
			cout<< "     : There is no restart file. Now, it will be initialized using the specified initial conditions." << endl;
#pragma omp parallel for num_threads(NT)
			for (int c = 0; c <= grid.NCM; c++)
			{
				grid.cells[c].data1D(indexV)	= IC_V(num_ic);
				grid.cells[c].data1D(indexQ)	= IC_Q(num_ic);
			}
		}
		else
		{
			int numVar;
			int ncm;
			bool is_read_OK = true;

			getline(restart_file, line);
			Common::DataRead::read_line(line);
			Common::DataRead::get_data_from_string<unsigned int>(line, "[ITERATION]", 0, iter);

			getline(restart_file, line);
			Common::DataRead::read_line(line);
			Common::DataRead::get_data_from_string<int>(line, "[NUMBER OF VARIABLES]:", 1, numVar);

			getline(restart_file, line);
			Common::DataRead::read_line(line);
			Common::DataRead::get_data_from_string<int>(line, "[NCM]:", 0, ncm);


			if (numVar != grid.cells[1].data1D.data[indexQ].data.size())
			{
				//throw Common::ExceptionDimensionMatch (FromHere(), "Number of variable in restart file does not match with the problem setting. Need to check the restart file");
				cout << "  [WARNING !!] Number of variable in restart file does not match with the problem setting. It will automatically use the IC of problem setting" << endl;
				cout <<	"				(Number of variable in restart file)   :  " << numVar << endl;
				cout <<	"				(Number of variable in problem setting):  " << grid.cells[1].data1D.data[indexQ].data.size() << endl;
				is_read_OK = false;
			}

			if (ncm != grid.NCM)
			{
				//throw Common::ExceptionDimensionMatch (FromHere(), "Number of cell in restart file does not match with the grid file. Need to check the restart file");
				cout << "  [WARNING !!] Number of cell in restart file does not match with the grid file. It will automatically use the IC of problem setting" << endl;
				cout <<	"				(Number of cell in a restart file)  :  " << ncm << endl;
				cout <<	"				(Number of cell in a grid file)     :  " << grid.NCM << endl;
				is_read_OK = false;
			}

			if (is_read_OK == true)
			{
				while (! restart_file.eof())
				{
					int c;
					restart_file >> c;

					for (int i = 0; i <= numVar-1; i++)	restart_file >> grid.cells[c].data1D(indexV).data[i];
				}

				restart_file.close();

#pragma omp parallel for num_threads(NT)
				for (int c = 1; c <= grid.NCM; c++)
				{
					CFD::VariableChange::V_to_Q(CFD_variabletype, CFD_NT, grid.cells[c].data1D(indexV), species_set, grid.ND, grid.cells[c].data1D(indexQ));
				}
			}
			else
			{
				restart_file.close();
#pragma omp parallel for num_threads(NT)
				for (int c = 0; c <= grid.NCM; c++)
				{
					grid.cells[c].data1D(indexV)	= IC_V(num_ic);
					grid.cells[c].data1D(indexQ)	= IC_Q(num_ic);
				}
			}
		}
	}
	else
	{
		cout<< ".....initializing with the specified ICs:" << endl;
#pragma omp parallel for num_threads(NT)
		for (int c = 0; c <= grid.NCM; c++)
		{
			grid.cells[c].data1D(indexV)	= IC_V(num_ic);
			grid.cells[c].data1D(indexQ)	= IC_Q(num_ic);
		}
	}



	if (problem_setup.is_viscous == true)
	{
		BCValuesWall.resize(problem_setup.NWCOND);

		for (int n = 0; n <= problem_setup.NWCOND-1; n++)
		{
			BCValuesWall(n)	= CFD::CFD_VariableSet::V(species_set, grid.ND, problem_setup.NER, problem_setup.NEV, problem_setup.NEE, problem_setup.is_viscous, problem_setup.is_axisymmetric);
		}

		string var_name_temp;

		for (int n = 0; n <= problem_setup.NWCOND-1; n++)
		{
			double rhoMix = 0.0;
			for (int s = 0; s <= species_set.NS-1; s++)
			{
				var_name_temp				= CFD::CFD_VariableSet::var_stringV(species_set.species[s].name, CFD::FluxCategory::Mass);
				BCValuesWall.data[n](var_name_temp) =  problem_setup.Ys[n][s];
				rhoMix	+= BCValuesWall.data[n](var_name_temp);
			}

			for (int s = 0; s <= species_set.NS-1; s++)
			{
				var_name_temp				= CFD::CFD_VariableSet::var_stringV(species_set.species[s].name, CFD::FluxCategory::Mass);
				BCValuesWall.data[n](var_name_temp) /= rhoMix;
			}


			for (int k = 0; k <= grid.ND-1; k++)
			{
				var_name_temp	= CFD::CFD_VariableSet::var_stringV(k, CFD::FluxCategory::Momentum);
				BCValuesWall.data[n](var_name_temp) =  0.0;
			}


			var_name_temp	=  CFD::CFD_VariableSet::var_stringV(0, CFD::FluxCategory::Energy);
			BCValuesWall.data[n](var_name_temp) =  problem_setup.Tw[n][0];

			if (problem_setup.NER != 0)
			{
				var_name_temp	=  CFD::CFD_VariableSet::var_stringV(CHEM::EnergyMode::E_ROT, CFD::FluxCategory::Energy);
				BCValuesWall.data[n](var_name_temp) =  problem_setup.Tw[n][1];
			}


			if (problem_setup.NEV != 0)
			{
				var_name_temp	=  CFD::CFD_VariableSet::var_stringV(CHEM::EnergyMode::E_VE, CFD::FluxCategory::Energy);
				BCValuesWall.data[n](var_name_temp) =  problem_setup.Tw[n][2];
			}


			if (problem_setup.NEE != 0)
			{
				var_name_temp	=  CFD::CFD_VariableSet::var_stringV(CHEM::EnergyMode::E_ELE, CFD::FluxCategory::Energy);
				BCValuesWall.data[n](var_name_temp) =  problem_setup.Tw[n][3];
			}
		}
	}


	DataTreatement(4, true, OP2A_DataCatogories::OP2A_CELL1D);
	DataTreatement(4, true, OP2A_DataCatogories::OP2A_CELL2D);
	DataTreatement(4, true, OP2A_DataCatogories::OP2A_FACE1D);
	DataTreatement(4, true, OP2A_DataCatogories::OP2A_FACE2D);
	check_elapsed_time("Initialization of Flow Conditions");
}

