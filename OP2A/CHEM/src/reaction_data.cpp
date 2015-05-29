/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Feb 11, 2015
 *      			Author: Minkwan Kim
 *
 * reaction_data.cpp
 * 			-  
 *  
 */




#include "../include/reaction.hpp"
#include "../include/OP2A_chemistry.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"




REACTION_DATA::REACTION_DATA()
{
	NR	= 0;
}

REACTION_DATA::~REACTION_DATA()
{
	if (data_entire.size() != 0)
		for (int i = 0; i <= data_entire.size()-1; i++)	data_entire[i]	= NULL;
}




REACTION_DATA_ver2::REACTION_DATA_ver2()
{
	NR	= 0;
	NS	= 0;
}

REACTION_DATA_ver2::~REACTION_DATA_ver2()
{

}


void REACTION_DATA_ver2::read_reaction_data(string file_name)
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

	NR	= 0;

	// READ REACTION DATA FILE
	if (reaction_file.is_open())
	{
		while (! reaction_file.eof())
		{
			getline(reaction_file, line);			data_read::read_line(line);
			if (line.compare(0, text_start.size(), text_start) == 0)	// START TO READ REACTION DATA
			{
				// Initialize variables
				reaction_k.push_back(REACTION_SINGLE());
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
					if(temp_string1.size() > 0)	reaction_k[NR].name	= temp_string1;

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

					// 3. READ PRODUCT INFO
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
					data_read::get_data_from_string<int>(line,	"TYPE:", -1, reaction_k[NR].type);


					// 5. FORWARD REACTION
					int temp;
					data_read::get_data_from_string<int>(line,	"FORWARD REACTION RATE METHOD:", 0, temp);
					if (line.compare(0, text_kf.size(), text_kf) == 0) sscanf(line.c_str(),"%s %s %s %lf %lf %lf %lf", temp_string2, temp_string3, temp_string4, &reaction_k[NR].kf_coeff[0],			 &reaction_k[NR].kf_coeff[1], 				&reaction_k[NR].kf_coeff[2], 			&reaction_k[NR].kf_coeff[3]);
					if (line.compare(0, text_Tf.size(), text_Tf) == 0) sscanf(line.c_str(),"%s %s %s %lf %lf %lf %lf", temp_string2, temp_string3, temp_string4, &reaction_k[NR].temperature_coeff_f[0], &reaction_k[NR].temperature_coeff_f[1], 	&reaction_k[NR].temperature_coeff_f[2], &reaction_k[NR].temperature_coeff_f[3]);

					// 6. BACKWARD REACTION
					data_read::get_data_from_string<int>(line,	"BACKWARD REACTION RATE METHOD:", 0, reaction_k[NR].method);
					if (line.compare(0, text_kb.size(), text_kb) == 0) sscanf(line.c_str(),"%s	%s %s %lf %lf %lf %lf", temp_string2, temp_string3, temp_string4, &reaction_k[NR].kb_coeff[0],			  &reaction_k[NR].kb_coeff[1], 				&reaction_k[NR].kb_coeff[2], 			&reaction_k[NR].kb_coeff[3]);
					if (line.compare(0, text_Tb.size(), text_Tb) == 0) sscanf(line.c_str(),"%s	%s %s %lf %lf %lf %lf", temp_string2, temp_string3, temp_string4, &reaction_k[NR].temperature_coeff_b[0], &reaction_k[NR].temperature_coeff_b[1], 	&reaction_k[NR].temperature_coeff_b[2], &reaction_k[NR].temperature_coeff_b[3]);


					// 7. Reference Temperature
					data_read::get_data_from_string<double>(line,	"REFERENCE TEMPERATURE:", 0.0, reaction_k[NR].Tref);


					// 8. EQUILIBRIUM CONSTANT DATA
					if (line.compare(0, text_Keq.size(), text_Keq) == 0)
					{
						data_read::get_data_from_string<int>(line,	"EQUILIBRIUM CONSTANTS:", 0, reaction_k[NR].Keq_data.model);
						getline(reaction_file, line);	data_read::read_line(line);

						int i	= 0;
						int nn = text_end.size();
						while (line.compare(0, text_end.size(), text_end) != 0)
						{
							reaction_k[NR].Keq_data.n.push_back(0.0);
							reaction_k[NR].Keq_data.A1.push_back(0.0);
							reaction_k[NR].Keq_data.A2.push_back(0.0);
							reaction_k[NR].Keq_data.A3.push_back(0.0);
							reaction_k[NR].Keq_data.A4.push_back(0.0);
							reaction_k[NR].Keq_data.A5.push_back(0.0);

							sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf", &reaction_k[NR].Keq_data.n[i],
																			&reaction_k[NR].Keq_data.A1[i],
																			&reaction_k[NR].Keq_data.A2[i],
																			&reaction_k[NR].Keq_data.A3[i],
																			&reaction_k[NR].Keq_data.A4[i],
																			&reaction_k[NR].Keq_data.A5[i]);
							i++;
							getline(reaction_file, line);	data_read::read_line(line);
						}
						reaction_k[NR].Keq_data.num_data	= i;
					}
				}


				// Update information
				reaction_k[NR].Reactant_coeff.resize(NS);
				reaction_k[NR].Product_coeff.resize(NS);
				for (int s = 0; s <= NS-1; s++)	reaction_k[NR].Reactant_coeff[s]	= 0.0;
				for (int s = 0; s <= NS-1; s++)	reaction_k[NR].Product_coeff[s]		= 0.0;


				for (int ss = 0; ss <= NS-1; ss++)
				{
					for (int s = 0; s <= num_reactant-1; s++)
					{
						if (react_name[s]	== species_data[ss].basic_data.name)
						{
							reaction_k[NR].Reactant_coeff[ss]	+= react_coeff[s];
						}
					}

					for (int s = 0; s <= num_product-1; s++)
					{
						if (prod_name[s]	== species_data[ss].basic_data.name)
						{
							reaction_k[NR].Product_coeff[ss]	+= prod_coeff[s];
						}
					}
				}

				NR++;
			}
		}
	}
	else
	{
		Error_message_type	error;
		error.location_primary_name 	= "reaction_data_read";
		error.location_secondary_name 	= "N/A";
		error.message					= "CANNOT FIND A REACTION DATA FILE. PLEASE CHECK THE FILE";
		error.print_message();
	}


	for (int k = 0; k <= NR-1; k++)
	{
		double sum_coeff_f = 0.0;
		double sum_coeff_b = 0.0;

		for (int i = 0; i <= 3; i++)
		{
			sum_coeff_f += reaction_k[k].temperature_coeff_f[i];
			sum_coeff_b += reaction_k[k].temperature_coeff_b[i];
		}

		for (int i = 0; i <= 3; i++)
		{
			reaction_k[k].temperature_coeff_f[i] /= sum_coeff_f;
			reaction_k[k].temperature_coeff_b[i] /= sum_coeff_b;
		}
	}


	Xf_all.resize(NR);
	Xf_molecule.resize(NR);
	Xf_electron.resize(NR);

	Xb_all.resize(NR);
	Xb_molecule.resize(NR);
	Xb_electron.resize(NR);


	for (int  k = 0; k <= NR-1; k++)
	{
		Xf_all[k].resize(NS);
		Xf_molecule[k].resize(NS);
		Xf_electron[k].resize(NS);

		Xb_all[k].resize(NS);
		Xb_molecule[k].resize(NS);
		Xb_electron[k].resize(NS);

		for (int s = 0; s <= NS-1; s++)
		{
			Xf_all[k][s]		= 0.0;
			Xf_molecule[k][s]	= 0.0;
			Xf_electron[k][s]	= 0.0;

			Xb_all[k][s]		= 0.0;
			Xb_molecule[k][s]	= 0.0;
			Xb_electron[k][s]	= 0.0;
		}
	}



	for (int k = 0; k <= NR-1; k++)
	{
		// FORWARD REACTION
		double sum_Xf_all		= 0.0;
		double sum_Xf_molecule	= 0.0;
		double sum_Xf_electron	= 0.0;

		for (int s = 0; s <= NS-1; s++)
		{
			switch(species_data[s].basic_data.type)
			{
			case MOLECULE:
				sum_Xf_all 		+= reaction_k[k].Reactant_coeff[s];
				sum_Xf_molecule += reaction_k[k].Reactant_coeff[s];

				Xf_all[k][s]		= reaction_k[k].Reactant_coeff[s];
				Xf_molecule[k][s]	= reaction_k[k].Reactant_coeff[s];
				break;

			case ATOM:
				sum_Xf_all 		+= reaction_k[k].Reactant_coeff[s];

				Xf_all[k][s]		= reaction_k[k].Reactant_coeff[s];
				break;
			case ELECTRON:
				sum_Xf_electron 		+= reaction_k[k].Reactant_coeff[s];

				Xf_electron[k][s]		= reaction_k[k].Reactant_coeff[s];
				break;
			}
		}

		for (int s = 0; s <= NS-1; s++)
		{
			if (sum_Xf_all > 0.0)		Xf_all[k][s]			/= sum_Xf_all;
			else						Xf_all[k][s]			= 0.0;

			if (sum_Xf_molecule > 0.0)	Xf_molecule[k][s]		/= sum_Xf_molecule;
			else						Xf_molecule[k][s]		= 0.0;

			if (sum_Xf_electron > 0.0)	Xf_electron[k][s]		/= sum_Xf_electron;
			else						Xf_electron[k][s]		= 0.0;
		}

		// BACKWARD REACTION
		double sum_Xb_all		= 0.0;
		double sum_Xb_molecule	= 0.0;
		double sum_Xb_electron	= 0.0;

		for (int s = 0; s <= NS-1; s++)
		{
			switch(species_data[s].basic_data.type)
			{
			case MOLECULE:
				sum_Xb_all 		+= reaction_k[k].Product_coeff[s];
				sum_Xb_molecule += reaction_k[k].Product_coeff[s];

				Xb_all[k][s]		= reaction_k[k].Product_coeff[s];
				Xb_molecule[k][s]	= reaction_k[k].Product_coeff[s];
				break;

			case ATOM:
				sum_Xb_all 		+= reaction_k[k].Product_coeff[s];

				Xb_all[k][s]	= reaction_k[k].Product_coeff[s];
				break;
			case ELECTRON:
				sum_Xb_electron 		+= reaction_k[k].Product_coeff[s];

				Xb_electron[k][s]		= reaction_k[k].Product_coeff[s];
				break;
			}
		}

		for (int s = 0; s <= NS-1; s++)
		{
			if (sum_Xb_all > 0.0)		Xb_all[k][s]			/= sum_Xb_all;
			else						Xb_all[k][s]			= 0.0;

			if (sum_Xb_molecule > 0.0)	Xb_molecule[k][s]		/= sum_Xb_molecule;
			else						Xb_molecule[k][s]		= 0.0;

			if (sum_Xb_electron > 0.0)	Xb_electron[k][s]		/= sum_Xb_electron;
			else						Xb_electron[k][s]		= 0.0;
		}
	}
}
