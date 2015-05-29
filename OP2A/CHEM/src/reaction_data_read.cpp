/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Feb 11, 2015
 *      			Author: Minkwan Kim
 *
 * reaction_data_read.cpp
 * 			-  
 *  
 */

#include "../include/reaction.hpp"
#include "../include/species_data_class.hpp"
#include "../include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"



void reaction_data_read(string file_name, REACTION_DATA &reactions, vector<SPECIES> &species_data)
{
	string		line;
	ifstream 	reaction_file;
	reaction_file.open(file_name.c_str());

	string text_start 		= "[REACTION START]";
	string text_end 		= "[REACTION END]";
	string text_reactant 	= "REACTANTS:";
	string text_product 	= "PRODUCTS:";
	string text_type 		= "TYPE:";
	string text_kf 			= "FORWARD REACTION COEFFICIENT:";
	string text_Tf 			= "FORWARD REACTION TEMPERATURE:";
	string text_kb 			= "BACKWARD REACTION COEFFICIENT:";
	string text_Tb 			= "BACKWARD REACTION TEMPERATURE:";
	string text_Keq			= "EQUILIBRIUM CONSTANTS:";


	string temp_string1;
	char   temp_string2[20];
	char   temp_string3[20];
	char   temp_string4[20];
	double temp_double1;

	int NR	= 0;
	int NS	= species_data.size();

	// READ REACTION DATA FILE
	if (reaction_file.is_open())
	{
		while (! reaction_file.eof())
		{
			getline(reaction_file, line);			data_read::read_line(line);
			if (line.compare(0, text_start.size(), text_start) == 0)	// START TO READ REACTION DATA
			{
				// Initialize variables
				reactions.data_entire.push_back(new REACTION());
				int	num_reactant	= 0;
				vector<int>			react_id;
				vector<double>		react_coeff;
				vector<string>		react_name;

				int num_product		= 0;
				vector<int>			prod_id;
				vector<double>		prod_coeff;
				vector<string>		prod_name;


				while (line.compare(0, text_end.size(), text_end) != 0)
				{
					getline(reaction_file, line);	data_read::read_line(line);

					// 1. READ NAME
					temp_string1	= data_read::get_data_from_string_string(line, "NAME:");
					if(temp_string1.size() > 0)	reactions.data_entire[NR]->name	= temp_string1;

					// 2. READ REACTANT INFO
					if (line.compare(0, text_reactant.size(), text_reactant) == 0)
					{
						while (line.compare(0, text_product.size(), text_product) != 0)
						{
							getline(reaction_file, line);	data_read::read_line(line);

							react_id.push_back(0);
							react_coeff.push_back(0.0);
							react_name.push_back("");

							temp_double1 = -1;
							sscanf(line.c_str(),"%lf	%s", &temp_double1, temp_string2);

							if (temp_double1 != -1)
							{
								react_coeff[num_reactant]	= temp_double1;
								react_name[num_reactant]	= temp_string2;
								num_reactant++;
							}
						}
					}

					// 3. READ REACTANT INFO
					if (line.compare(0, text_product.size(), text_product) == 0)
					{
						while (line.compare(0, text_type.size(), text_type) != 0)
						{
							getline(reaction_file, line);	data_read::read_line(line);

							prod_id.push_back(0);
							prod_coeff.push_back(0.0);
							prod_name.push_back("");

							temp_double1 = -1;
							sscanf(line.c_str(),"%lf	%s", &temp_double1, temp_string2);

							if (temp_double1 != -1)
							{
								prod_coeff[num_product]	= temp_double1;
								prod_name[num_product]	= temp_string2;
								num_product++;
							}
						}
					}


					// 4. READ REACTION TYPE
					data_read::get_data_from_string<int>(line,	"TYPE:", -1, reactions.data_entire[NR]->type);


					// 5. FORWARD REACTION
					data_read::get_data_from_string<int>(line,	"FORWARD REACTION RATE METHOD:", 0, reactions.data_entire[NR]->kf.method);
					if (line.compare(0, text_kf.size(), text_kf) == 0) sscanf(line.c_str(),"%s %s %s %lf %lf %lf %lf", temp_string2, temp_string3, temp_string4, &reactions.data_entire[NR]->kf.rate_coeff[0],	&reactions.data_entire[NR]->kf.rate_coeff[1], 	&reactions.data_entire[NR]->kf.rate_coeff[2], 	&reactions.data_entire[NR]->kf.rate_coeff[3]);
					if (line.compare(0, text_Tf.size(), text_kf) == 0) sscanf(line.c_str(),"%s %s %s %lf %lf %lf %lf", temp_string2, temp_string3, temp_string4, &reactions.data_entire[NR]->kf.temperature[0], &reactions.data_entire[NR]->kf.temperature[1], 	&reactions.data_entire[NR]->kf.temperature[2], 	&reactions.data_entire[NR]->kf.temperature[3]);


					// 6. BACKWARD REACTION
					data_read::get_data_from_string<int>(line,	"BACKWARD REACTION RATE METHOD:", 0, reactions.data_entire[NR]->kb.method);
					if (line.compare(0, text_kb.size(), text_kb) == 0)	sscanf(line.c_str(),"%s	%s %s %lf %lf %lf %lf", temp_string2, temp_string3, temp_string4, &reactions.data_entire[NR]->kb.rate_coeff[0],	&reactions.data_entire[NR]->kb.rate_coeff[1], 	&reactions.data_entire[NR]->kb.rate_coeff[2], 	&reactions.data_entire[NR]->kb.rate_coeff[3]);
					if (line.compare(0, text_Tb.size(), text_kb) == 0)	sscanf(line.c_str(),"%s	%s %s %lf %lf %lf %lf", temp_string2, temp_string3, temp_string4, &reactions.data_entire[NR]->kb.temperature[0],&reactions.data_entire[NR]->kb.temperature[1], 	&reactions.data_entire[NR]->kb.temperature[2], 	&reactions.data_entire[NR]->kb.temperature[3]);


					// 7. Reference Temperature
					data_read::get_data_from_string<double>(line,	"REFERENCE TEMPERATURE:", 0.0, reactions.data_entire[NR]->Tref);


					// 8. EQUILIBRIUM CONSTANT DATA
					if (line.compare(0, text_Keq.size(), text_Keq) == 0)
					{
						data_read::get_data_from_string<int>(line,	"EQUILIBRIUM CONSTANTS:", 0, reactions.data_entire[NR]->Keq.model);
						getline(reaction_file, line);	data_read::read_line(line);

						int i	= 0;
						int nn = text_end.size();
						while (line.compare(0, text_end.size(), text_end) != 0)
						{
							reactions.data_entire[NR]->Keq.n.push_back(0.0);
							reactions.data_entire[NR]->Keq.A1.push_back(0.0);
							reactions.data_entire[NR]->Keq.A2.push_back(0.0);
							reactions.data_entire[NR]->Keq.A3.push_back(0.0);
							reactions.data_entire[NR]->Keq.A4.push_back(0.0);
							reactions.data_entire[NR]->Keq.A5.push_back(0.0);

							sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf", &reactions.data_entire[NR]->Keq.n[i],
																			&reactions.data_entire[NR]->Keq.A1[i],
																			&reactions.data_entire[NR]->Keq.A2[i],
																			&reactions.data_entire[NR]->Keq.A3[i],
																			&reactions.data_entire[NR]->Keq.A4[i],
																			&reactions.data_entire[NR]->Keq.A5[i]);
							i++;
							getline(reaction_file, line);	data_read::read_line(line);
						}
						reactions.data_entire[NR]->Keq.num_data	= i;
					}
				}


				// Update reactant ID information
				reactions.data_entire[NR]->Reactant_num	= num_reactant;
				reactions.data_entire[NR]->Reactant_coeff.resize(num_reactant);
				reactions.data_entire[NR]->Reactant_id.resize(num_reactant);

				for (int s = 0; s <= num_reactant-1; s++)
				{
					reactions.data_entire[NR]->Reactant_coeff[s]	= react_coeff[s];
					for (int ss = 0; ss <= NS-1; ss++)
					{
						if (react_name[s]	== species_data[ss].basic_data.name)
						{
							reactions.data_entire[NR]->Reactant_id[s]	= ss;
						}
					}
				}

				// Update Product ID
				reactions.data_entire[NR]->Product_num	= num_product;
				reactions.data_entire[NR]->Product_coeff.resize(num_reactant);
				reactions.data_entire[NR]->Product_id.resize(num_reactant);

				for (int s = 0; s <= num_product-1; s++)
				{
					reactions.data_entire[NR]->Product_coeff[s]	= prod_coeff[s];
					for (int ss = 0; ss <= NS-1; ss++)
					{
						if (prod_name[s]	== species_data[ss].basic_data.name)
						{
							reactions.data_entire[NR]->Product_id[s]	= ss;
						}
					}
				}

				NR++;

			}
		}

		reactions.NR	= NR;
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "reaction_data_read";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A REACTION DATA FILE. PLEASE CHECK THE FILE";
		error.print_message();
	}

}









