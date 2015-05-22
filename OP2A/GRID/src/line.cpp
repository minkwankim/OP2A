/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Mar 16, 2015
 *      			Author: Minkwan Kim
 *
 * line.cpp
 * 			-  
 *  
 */


#include <stddef.h>
#include "../include/grid_class.hpp"
#include "../include/cell_class.hpp"
#include "../include/constants_grid.hpp"
#include "../../UTIL/include/OP2A_utilities.hpp"


bool determine_line_starting_point(unsigned int BC)
{
	bool start;
	if (BC == BC_WALL ||	BC == BC_ANODE ||	BC == BC_CATHODE ||	BC == BC_DIELECTRIC_WALL)	start = true;
	else																						start = false;

	return(start);
}


void line_finder(GRID_CLASS *grid, vector <vector <int> > &lines, vector <vector <int> > &lines_bd, vector <int> &cell_line_ID, int &line_num)
{
	lines_bd	= vector_2D(grid->NCM+1, 2, 0);
	cell_line_ID.resize(grid->NCM+1);
	for (int i = 0; i <= grid->NCM; i++)	cell_line_ID[i] = 0;



	/*
	 * 1. Find number of walls
	 * 		- Assign Face IDs on the starting point of lines
	 */
	int n_wall_faces = 0;
	int	wall_face_ptr;
	int cell_ptr;
	vector <int> wall;	// Face IDs of starting points
	wall.push_back(-1);

	for (int c = 0; c <= grid->NGM-1; c++)
	{
		if (determine_line_starting_point(grid->cells_ghost.data_ptr[c]->BC) == true)
		{
			n_wall_faces++;
			wall.push_back(grid->cells_ghost.data_ptr[c]->face[0]);
		}
	}

	// ERROR CHECK
	if (n_wall_faces > MAX_NUM_LINE)
	{
		Error_message_type	error;
		error.message		= "Number of lines exceed the MAX_NUM_LINE. Need to increase NUM_LINE in constants_grid->hpp";
		error.module_name	= "GRID: line.cpp";
		error.print_message();
	}



	/*
	 * 2. Initialize variables
	 * 		- lines_new [line_number][0] ( 0: open / 1:close)
	 * 								 [position in line] (Cell IDs (not pointer))
	 *
	 * 		- line_optiions[line_number][0] (Total number of available options
	 * 									[ ] (list of options: cell IDs)
	 *
	 * 		- cell_line_ID[cell_ID]	(line ID number of each cells)
	 */
	vector <vector <int> > lines_new		= vector_2D<int> (n_wall_faces+1, MAX_LINE_LENGTH+1, 0);
	vector <vector <int> > line_options		= vector_2D<int> (n_wall_faces+1, C_MAX_FACES+1, 0);
	vector <double>	direction(grid->DIM, 0.0);
	vector <double>	trial(grid->DIM, 0.0);




	/*
	 * 3. Assign open lines to starting points
	 * 		- Allocate starting cells on all lines
	 */
	int open_line = 0;		// number of open lines
	int f, cl;
	int cc, cc_ptr, f_ptr;

	int pos = 1;
	for (int c = 0; c <= grid->NGM-1; c++)
	{
		if (determine_line_starting_point(grid->cells_ghost.data_ptr[c]->BC) == true)
		{
			f	= grid->faces.whereis[grid->cells_ghost.data_ptr[c]->face[0]];
			cl	= grid->faces.data_ptr[f]->cl[0];

			open_line++;
			cell_line_ID[cl]			= open_line;
			lines_new[open_line][pos]	= cl;
			lines_new[open_line][0] 	= 1;	// Open line
		}
	}



	/*
	 * 4. Build lines
			- line candidates [cell_ID][possible options] (line number candidate)
	 */
	int total_options;
	int candidate;

	while (open_line != 0)
	{
		pos++;
		if (pos > MAX_LINE_LENGTH)
		{
			Error_message_type	error;
			error.message		= "Length of lines exceed the MAX_LINE_LENGTH. Need to increase NUM_LINE_LENGTH in constants_grid->hpp";
			error.module_name	= "GRID: line.cpp";
			error.print_message();
		}

		// Reset counters
		open_line		= 0;
		total_options 	= 0;
		vector <vector <int> > line_candidates	= vector_2D<int> (grid->NCM+1, C_MAX_FACES+1, 0);

		// Find all possible options to glow line
		for (int s = 1; s <= n_wall_faces; s++)
		{
			wall_face_ptr	= grid->faces.whereis[wall[s]];

			if (lines_new[s][0]	== 1)	// Check whether line is opened
			{
				lines_new[s][0]	= 0;	// Temporary close line
				cell_ptr		= grid->cells.whereis[lines_new[s][pos-1]];	// Get cell of previous position in line

				for (int k = 0; k <= grid->DIM-1; k++)	direction[k]	= -grid->faces.data_ptr[wall_face_ptr]->n[0][k];

				int n_options	= 0;
				for (int f = 0; f <= grid->cells.data_ptr[cell_ptr]->NF-1; f++)
				{
					int ff	= grid->cells.data_ptr[cell_ptr]->face[f];
					f_ptr	= grid->faces.whereis[ff];
					cc		= grid->faces.data_ptr[f_ptr]->cl[0] + grid->faces.data_ptr[f_ptr]->cr[0] - lines_new[s][pos-1];

					if (cc > 0 && cell_line_ID[cc] <= 0)
					{
						// Calculate dot product
						cc_ptr = grid->cells.whereis[cc];
						for (int k = 0; k <= grid->DIM-1; k++)	trial[k] = grid->cells.data_ptr[cc_ptr]->x[k]	- grid->cells.data_ptr[cell_ptr]->x[k];

						double length = 0.0;
						for (int k = 0; k <= grid->DIM-1; k++)	length		+= trial[k]*trial[k];

						length	= sqrt(length);
						for (int k = 0; k <= grid->DIM-1; k++)	trial[k]	/= length;

						double dot = 0.0;
						for (int k = 0; k <= grid->DIM-1; k++)	dot += trial[k]*direction[k];


						if (dot > 0.0)
						{
							n_options++;
							line_options[s][0]			= n_options;
							line_options[s][n_options]	= cc;

							candidate	= line_candidates[cc][0] + 1;
							line_candidates[cc][0]			= candidate;
							line_candidates[cc][candidate] 	= s;
						}
					}
				}
			}
		}

		// Count total number of options and open lines
		for (int c = 1; c <= grid->NCM; c++)	if (line_candidates[c][0] > 0) total_options++;
		for (int s = 1; s <= n_wall_faces; s++)	if (line_options[s][0] > 0)	open_line++;


		while (total_options > 0 && open_line > 0)
		{
			int focus = C_MAX_FACES;

			// Find minimum number of options and start from it
			for (int s = 1; s <= n_wall_faces; s++)	if (line_options[s][0] > 0)	focus	= fmin(line_options[s][0], focus);

			for (int s = 1; s <= n_wall_faces; s++)
			{
				if(line_options[s][0] == focus)
				{
					cc = line_options[s][1];

					if (cell_line_ID[cc] <= 0)	// Line is not assigned yet
					{
						cell_line_ID[cc] 	= s;
						lines_new[s][0] 	= 1;
						lines_new[s][pos]	= cc;
						total_options--;

						// Initialize selected line option
						for (int j = 1; j <= line_options[s][0]; j++)	line_options[s][j]	= 0;
						line_options[s][0]	= 0;

						for (int j = 1; j <= line_candidates[cc][0]; j++)
						{
							remove_element_in_vector<int>(cc, line_options[line_candidates[cc][j]]);	// Remove selectred cell from the list

							// Reset candidate of the selected cell
							line_candidates[cc][0]	= 0;
							line_candidates[cc][j]	= 0;
						}
					}
				}
			}

			// Find number of opened lines
			open_line = 0;
			for (int s = 1; s <= n_wall_faces; s++)	if (line_options[s][0] > 0) open_line++;
		}

		open_line = 0;
		for (int s = 1; s <= n_wall_faces; s++)	if (lines_new[s][0] > 0)	open_line++;
	}


	/*
	 * 5. Make it a single array
	 * 		- lines[line number][0] (start position in line data base)
	 * 						   [1] (end position in line data base)
	 *
	 * 		- line_bd[position number][0] = cell_ID
	 * 								  [1] = face_ID
	 */
	lines	= vector_2D<int>(n_wall_faces+2, 2, 0);
	int num_lines 	= 0;
	int	line_length	= 0;
	int k;
	pos	= 0;


	for (int s = 1; s <= n_wall_faces; s++)
	{
		// 5.1 Get length of line
		line_length = 0;
		while(lines_new[s][line_length+1] != 0)	line_length++;

		// 5.2 Assign the informatio into a single array
		if (line_length > 1)
		{
			num_lines++;

			k = 1;						// Set starting position in line
			int c = lines_new[s][k];	// Get cell ID

			// Start building array
			if (c != 0)
			{
				lines[num_lines][0]	= pos;			// Starting position in data array

				// Set starting data of line
				lines_bd[pos][0]	= 	c;			// Cell ID
				lines_bd[pos][1]	=	-1;			// Face ID

				k++;
				pos++;
				cc	= lines_new[s][k];

				// Get data of next cell/face on the line
				//	- c  = previous cell
				//	- cc = next cell
				while (cc != 0)
				{
					int found_face 	= 0;
					int face_current;
					int	face_next;

					for (int f = 0; f <= grid->cells.data_ptr[grid->cells.whereis[c]]->NF-1; f++)
					{
						face_current	= grid->cells.data_ptr[grid->cells.whereis[c]]->face[f];
						if ((grid->faces.data_ptr[grid->faces.whereis[face_current]]->cl[0] + grid->faces.data_ptr[grid->faces.whereis[face_current]]->cr[0] - c) == cc)
						{
							face_next = face_current;
							found_face++;
						}
					}

					// CHECK ERROR
					if (found_face != 1)
					{
						Error_message_type	error;
						error.message		= "Problem in finding face AT constants_grid->hpp";
						error.module_name	= "GRID: line.cpp";
						error.print_message();
					}

					// Update end position and data set
					lines[num_lines][1]	= pos;
					lines_bd[pos][0]	= cc;
					lines_bd[pos][1]	= face_next;

					pos++;
					k++;

					c	= cc;
					cc 	= lines_new[s][k];
				}
			}
			else
			{
				Error_message_type	error;
				error.message		= "Line does not start correctly!";
				error.module_name	= "GRID: line.cpp";
				error.print_message();
			}
		}
		else
		{
			cout << "Line #" << s << "has " << line_length << " cells, so it will be deleted" << endl;
			for (int j = 1; j <= line_length; j++)
			{
				int c	= lines_new[s][j];
				cell_line_ID[c] = -1;
			}
		}
	}

	// 5.3. Account cells don't belong to any lines
	lines[num_lines+1][0]	= 	pos;	// START
	lines[num_lines+1][1]	=	0;		// END

	for (int c = 1; c <= grid->NCM; c++)
	{
		if (cell_line_ID[c]	== -1)
		{
			lines_bd[pos][0]		= c;
			lines_bd[pos][1]		= -1;

			lines[num_lines+1][1]	= pos;
			pos++;
		}
	}
}



