/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Feb 6, 2015
 *      			Author: Minkwan Kim
 *
 * problem_WALLCOND.cpp
 * 			-  
 *  
 */




#include "../../UTIL/include/OP2A_utilities.hpp"
#include "../include/OP2A_problem_wallcond.hpp"
#include "../include/OP2A_setup_constants.hpp"


PROBLEM_WALLCOND::PROBLEM_WALLCOND()
{
	NWCOND	= 0;
}

PROBLEM_WALLCOND::~PROBLEM_WALLCOND()
{

}



void PROBLEM_WALLCOND::read_WALLCOND(string filename)
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
	NWCOND	= 0;
	if (prob_file.is_open())
	{
		while (! prob_file.eof())
		{
			getline(prob_file, line);
			data_read::read_line(line);

			temp_text = "#WALL COND NUM[START]: ";
			index_temp = -1;
			data_read::get_data_from_string<int>(line, temp_text, -1, index_temp);

			if (index_temp >= 0)
			{
				Ys.push_back(vector<double>(MAX_NS_PROBLEM, 0.0));
				Tw.push_back(vector<double>(4, 0.0));
				emissivity_Tref.push_back(0.0);
				emissivity_below.push_back(vector<double>(6, 0.0));
				emissivity_above.push_back(vector<double>(6, 0.0));


				temp_text1 = "#INFLOW COND NUM[END]";
				while (line.compare(0, temp_text1.size(), temp_text1) != 0)
				{
					getline(prob_file, line);
					data_read::read_line(line);

					temp_text 	= "Ys[";
					flag1	  	= false;
					s_temp		= -1;
					data_read::get_data_from_string<int>(line, temp_text, -1, s_temp);

					// Mass fraction
					if (s_temp != -1)
					{
						temp_text2 = temp_text;
						temp_text3 = convertInt(s_temp);
						temp_text2 += temp_text3;
						temp_text2 += "]:";
						data_read::get_data_from_string<double>(line, temp_text2, 1.0e-15, temp1);
						Ys[NWCOND][s_temp]	= temp1;
					}

					// Temperature
					temp1 = -1.0;
					data_read::get_data_from_string<double>(line, "Tw:", 0.0, temp1);
					if (temp1 != -1.0)	Tw[NWCOND][0]	= temp1;

					temp1 = -1.0;
					data_read::get_data_from_string<double>(line, "Twr:", 0.0, temp1);
					if (temp1 != -1.0)	Tw[NWCOND][1]	= temp1;

					temp1 = -1.0;
					data_read::get_data_from_string<double>(line, "Twv:", 0.0, temp1);
					if (temp1 != -1.0)	Tw[NWCOND][2]	= temp1;

					temp1 = -1.0;
					data_read::get_data_from_string<double>(line, "Twe:", 0.0, temp1);
					if (temp1 != -1.0)	Tw[NWCOND][3]	= temp1;

					data_read::get_data_from_string<double>(line, "EMISSIVITY COEFFICIENT REF TEMPERATURE: ", 0.0, emissivity_Tref[NWCOND]);
					data_read::get_data_from_string<double>(line, "EMISSIVITY COEFFICIENT BELOW REF TEMPERATURE A:", 0.0, emissivity_below[NWCOND][0]);
					data_read::get_data_from_string<double>(line, "EMISSIVITY COEFFICIENT BELOW REF TEMPERATURE B:", 0.0, emissivity_below[NWCOND][1]);
					data_read::get_data_from_string<double>(line, "EMISSIVITY COEFFICIENT BELOW REF TEMPERATURE C:", 0.0, emissivity_below[NWCOND][2]);
					data_read::get_data_from_string<double>(line, "EMISSIVITY COEFFICIENT BELOW REF TEMPERATURE D:", 0.0, emissivity_below[NWCOND][3]);
					data_read::get_data_from_string<double>(line, "EMISSIVITY COEFFICIENT BELOW REF TEMPERATURE E:", 0.0, emissivity_below[NWCOND][4]);
					data_read::get_data_from_string<double>(line, "EMISSIVITY COEFFICIENT BELOW REF TEMPERATURE F:", 0.0, emissivity_below[NWCOND][5]);

					data_read::get_data_from_string<double>(line, "EMISSIVITY COEFFICIENT ABOVE REF TEMPERATURE A:", 0.0, emissivity_above[NWCOND][0]);
					data_read::get_data_from_string<double>(line, "EMISSIVITY COEFFICIENT ABOVE REF TEMPERATURE B:", 0.0, emissivity_above[NWCOND][1]);
					data_read::get_data_from_string<double>(line, "EMISSIVITY COEFFICIENT ABOVE REF TEMPERATURE C:", 0.0, emissivity_above[NWCOND][2]);
					data_read::get_data_from_string<double>(line, "EMISSIVITY COEFFICIENT ABOVE REF TEMPERATURE D:", 0.0, emissivity_above[NWCOND][3]);
					data_read::get_data_from_string<double>(line, "EMISSIVITY COEFFICIENT ABOVE REF TEMPERATURE E:", 0.0, emissivity_above[NWCOND][4]);
					data_read::get_data_from_string<double>(line, "EMISSIVITY COEFFICIENT ABOVE REF TEMPERATURE F:", 0.0, emissivity_above[NWCOND][5]);
				}
				NWCOND++;
			}
		}
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "OP2A_problem_WALLCOND";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A PROBLEM SETUP DATA FILE. PLEASE CHECK THE FILE";
		error.print_message();
	}

	prob_file.close();
	prob_file.clear();
}
