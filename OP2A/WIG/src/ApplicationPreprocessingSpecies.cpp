/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 12, 2015
 *      			Author: Minkwan Kim
 *
 * ApplicationPreprocessingSpecies.cpp
 * 			-  
 *  
 */




#include "Common/include/CodeLocation.hpp"
#include "Common/include/Exception_NPExceed.hpp"
#include "Common/include/StringOps.hpp"

#include "../include/OP2A_Application.hpp"
#include "CFD/include/VariableChange.hpp"



void ApplicationOP2A::preprocessing_species()
{
	 species_set.read_SpeciesSet(problem_setup.species_file, problem_setup.NS);
	 species_set.showInfo();
	 check_elapsed_time("Reading Species data set");

	 create_sampleDataCFD();
	 CFD_variabletype = CFD::VariableChange::VariableType(problem_setup.NER, problem_setup.NEV, problem_setup.NEE);
}
