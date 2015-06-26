/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 25, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_Application.hpp
 * 			-  
 *  
 */
#ifndef OP2A_APPLICATION_HPP_
#define OP2A_APPLICATION_HPP_


#include <iostream>
#include <fstream>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#include "Common/include/OP2A.hpp"
#include "Common/include/OP2A_time.hpp"
#include "Common/include/Version.hpp"

#include "DATA/include/DataStorage.hpp"
#include "DATA/include/DataStorage2D.hpp"
#include "DATA/include/DataStorageVector.hpp"

#include "GRID/include/Grid.hpp"
#include "GRID/include/PrintResult.hpp"

#include "CHEM/include/SpeciesSet.hpp"

#include "../include/OP2A_Problem.hpp"

using namespace std;
using namespace OP2A;
using namespace OP2A::Common;


#define	OP2A_VERSION_MAIN	1
#define	OP2A_VERSION_SUB	1
#define OP2A_MAX_N_TASK		100



enum	OP2A_ParallelMethod
{
	OP2A_NONE	= 0,
	OP2A_MPI	= 1,
	OP2A_OPENMP	= 2,
	OP2A_HYBRID	= 3
};

enum	OP2A_ParallelProcessor
{
	OP2A_CPU			= 0,
	OP2A_GPU			= 1,
	OP2A_MIC			= 2,
	OP2A_HYBRID_CPUGPU	= 3,
	OP2A_HYBRID_CPUMIC	= 4,
	OP2A_HYBRID_GPUMIC	= 5
};


enum	OP2A_DataCatogories
{
	OP2A_CELL1D			= 0,
	OP2A_CELL2D			= 1,
	OP2A_FACE1D			= 2,
	OP2A_FACE2D			= 3
};


class ApplicationOP2A
{
public:
	OP2A_ParallelMethod		parallel_method;
	OP2A_ParallelProcessor	parallel_processor;

	unsigned int	iter;
	int			P;
	int			NT;
	int			NP;
	bool		FLAG_TERMINATION;

	double 		t0;
	double 		t_simulation;
	double 		dt;

	const string		setup_filename;

	CPUTime time_running;						// Simulation time


	unsigned int CFD_NT;
	unsigned int CFD_variabletype;

	bool				use_extended_Stencil;
	OP2A_PROBLEM 		problem_setup;
	CHEM::SpeciesSet	species_set;
	GRID::Grid			grid;



	/*
	 * Constructor
	 * @author	Minkwan Kim
	 * @version	1.0	25/5/2015
	 */
	explicit ApplicationOP2A(OP2A_ParallelMethod app_parallel_method, OP2A_ParallelProcessor app_parallel_processor, int app_NT, string app_setup_filename)
						: parallel_method(app_parallel_method), parallel_processor(app_parallel_processor), NT(app_NT), setup_filename(app_setup_filename)
	{
		iter	= 0;
		P	= 0;
		NP	= 1;
		CFD_NT	= NT;
		CFD_variabletype = 0;
		FLAG_TERMINATION = false;

		t0	= 0.0;
		t_simulation = 0.0;
		dt = 0.0;

		m_t	= time(0);
		m_now = localtime(& m_t);

		use_extended_Stencil = false;

		indexQ		= 0;
		indexV		= 1;
		indexW		= 2;
		indexMIX	= 3;
		indexXs		= 4;
		indexYs		= 5;


		termination = false;
		RHS_n		= 0;
		RHS_max		= 0.0;
		RHS_2		= 0.0;
		CFLNumber	= 0.001;
	}

	~ApplicationOP2A()
	{

	}




private:
	string	m_module_name;

	time_t	m_t;
	struct	tm * m_now;


protected:
	Data::DataStorageVector<Data::DataStorage>		cell_data1D_template;
	Data::DataStorageVector<Data::DataStorage2D>	cell_data2D_template;
	Data::DataStorageVector<Data::DataStorage>		face_data1D_template;
	Data::DataStorageVector<Data::DataStorage2D>	face_data2D_template;

	Data::DataStorageVector<Data::DataStorage>		IC_Q;
	Data::DataStorageVector<Data::DataStorage>		IC_V;
	Data::DataStorageVector<Data::DataStorage>		BCValuesWall;

	int indexQ;
	int indexV;
	int indexW;
	int indexMIX;
	int indexXs;
	int indexYs;


public:
	bool 	termination;
	int 	RHS_n;
	double 	RHS_max;
	double 	RHS_2;
	double	CFLNumber;


public:
	void preparation(int argc, char *argv[], string modulename);
	void check_elapsed_time(string workname);
	void show_starting_task(string workname);

	void create_sampleDataCFD();

	void preprocessing_species();
	void preprocessing_grid();

	void print_result(const string& i_variablename, int type);
	void print_restartCFD(const string& i_variablename);

	void InitializeData(unsigned int num_ic, bool use_restart_file);
	void DataTreatement(int typeCase, bool is_initialize, int type);

	void ApplyBCInviscidNormal();
	void ApplyBCViscousNormal();

	void PrintConvergences(bool firststart);
	void CalcualteCFL();
	void Calcualtedt();


protected:
	void CalculateIC();
	void Cell1DDataTreatement(int typeCase, bool is_initialize);
	void Cell2DDataTreatement(int typeCase, bool is_initialize);
	void Face1DDataTreatement(int typeCase, bool is_initialize);
	void Face2DDataTreatement(int typeCase, bool is_initialize);




};




#endif /* OP2A_APPLICATION_HPP_ */
