/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jan 14, 2015
 *      			Author: Minkwan Kim
 *
 * problem_IC.cpp
 * 			-  
 *  
 */


#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../include/OP2A_problem_IC.hpp"
#include "../include/OP2A_setup_constants.hpp"


PROBLEM_IC_BASIC::PROBLEM_IC_BASIC()
{
	INITIALIZE_METHOD	= 0;						/* METHOD FOR INITIALIZATION => 0: USING INFLOW CONDITIONS / 1: USING AMBIENT CONDITIONS */
	NIC					= 0;
}

PROBLEM_IC_BASIC::~PROBLEM_IC_BASIC()
{

}




void PROBLEM_IC_BASIC::read_IC(string filename)
{



	int		temp_int;
	int 		index_temp;
	int 		s_temp, flag1;
	double 		temp1, temp2;
	string		temp_text;
	string 		temp_text1;
	string 		temp_text2;
	string 		temp_text3;
	string		line;
	ifstream	prob_file;

	prob_file.open(filename.c_str());

	// READ PROBLEM SETUP DATA FILE
	NIC		= 0;
	if (prob_file.is_open())
	{
		while (! prob_file.eof())
		{
			getline(prob_file, line);
			data_read::read_line(line);

			temp_int	= -1;
			data_read::get_data_from_string<int>(line, "INITIALIZE METHOD:", 0, temp_int);
			if (temp_int != -1)	INITIALIZE_METHOD	= temp_int;

			temp_text = "#INFLOW COND NUM[START]: ";
			index_temp = -1;
			data_read::get_data_from_string<int>(line, temp_text, -1, index_temp);

			if (index_temp >= 0)
			{
				rho_s.push_back(vector<double>(MAX_NS_PROBLEM, 0.0));
				v.push_back(vector<double>(3, 0.0));
				T.push_back(0.0);
				Tr.push_back(0.0);
				Tv.push_back(0.0);
				Te.push_back(0.0);

				/*
				T_wall.push_back(0.0);
				Tr_wall.push_back(0.0);
				Tv_wall.push_back(0.0);
				Te_wall.push_back(0.0);
				*/

				temp_text1 = "#INFLOW COND NUM[END]";
				while (line.compare(0, temp_text1.size(), temp_text1) != 0)
				{
					getline(prob_file, line);
					data_read::read_line(line);

					temp_text 	= "RHO[";
					flag1	  	= false;
					s_temp		= -1;
					data_read::get_data_from_string<int>(line, temp_text, -1, s_temp);

					// SPECIES Density
					if (s_temp != -1)
					{
						temp_text2 = temp_text;
						temp_text3 = convertInt(s_temp);
						temp_text2 += temp_text3;
						temp_text2 += "]:";
						data_read::get_data_from_string<double>(line, temp_text2, 1.0e-15, temp1);
						rho_s[NIC][s_temp]	= temp1;
					}

					// Velocity
					data_read::get_data_from_string<double>(line, "Vx: ", 0.0, v[NIC][0]);
					data_read::get_data_from_string<double>(line, "Vy: ", 0.0, v[NIC][1]);
					data_read::get_data_from_string<double>(line, "Vz: ", 0.0, v[NIC][2]);


					// Temperature
					temp1 = -1.0;
					data_read::get_data_from_string<double>(line, "T:", 0.0, temp1);
					if (temp1 != -1.0)	T[NIC]	= temp1;

					temp1 = -1.0;
					data_read::get_data_from_string<double>(line, "Tr:", 0.0, temp1);
					if (temp1 != -1.0)	Tr[NIC]	= temp1;

					temp1 = -1.0;
					data_read::get_data_from_string<double>(line, "Tv:", 0.0, temp1);
					if (temp1 != -1.0)	Tv[NIC]	= temp1;

					temp1 = -1.0;
					data_read::get_data_from_string<double>(line, "Te:", 0.0, temp1);
					if (temp1 != -1.0)	Te[NIC]	= temp1;
				}
				NIC++;
			}
		}
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "OP2A_problem_IC";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A PROBLEM SETUP DATA FILE. PLEASE CHECK THE FILE";
		error.print_message();
	}

	rho.resize(NIC);

	prob_file.close();
	prob_file.clear();


	whereis.resize(NIC);
	prob_file.open(filename.c_str());

	// READ PROBLEM SETUP DATA FILE
	NIC		= 0;
	if (prob_file.is_open())
	{
		while (! prob_file.eof())
		{
			getline(prob_file, line);
			data_read::read_line(line);

			temp_text = "#INFLOW COND NUM[START]: ";
			index_temp = -1;
			data_read::get_data_from_string<int>(line, temp_text, -1, index_temp);

			if (index_temp >= 0)
			{
				whereis[index_temp]	= NIC;
				NIC++;
			}
		}
	}
	prob_file.close();
}





PROBLEM_IC_CFD::PROBLEM_IC_CFD()
{

}


PROBLEM_IC_CFD::~PROBLEM_IC_CFD()
{

}


void PROBLEM_IC_CFD::assign_size()
{
	if (NIC > 0)
	{
		E.resize(NIC);
		Er.resize(NIC);
		Ev.resize(NIC);
		Ee.resize(NIC);
	}
}
