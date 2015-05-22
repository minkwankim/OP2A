/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Apr 14, 2015
 *      			Author: Minkwan Kim
 *
 * OP2A_read_write.hpp
 * 			-  
 *  
 */

#ifndef OP2A_READ_WRITE_HPP_
#define OP2A_READ_WRITE_HPP_



#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>


using namespace std;

namespace data_read
{
	string get_data_from_string_string(string line, string TEXT);

	string remove_space_front(string line);
	string remove_space_end(string line);

	void remove_comments(string& line, string comment_s, string comment_e);
	void remove_comments(string& line, string comment_s);

	void read_line(string& line);

	template<class T>
	void get_data_from_string(string line, string TEXT, T default_value, T& data)
	{
		int num;
		T result;

		// REMOVE EMPTY SPACE AT FRONT AND END
		line = remove_space_front(line);
		line = remove_space_end(line);

		// REMOVE EMPTY SPACE AT FRONT AND END
		TEXT = remove_space_front(TEXT);
		TEXT = remove_space_end(TEXT);

		num = TEXT.size();
		if(line.compare(0, num, TEXT) == 0)
		{
			line.erase(0, num);
			line = remove_space_front(line);
			istringstream(line) >> result;

			if (line == "") result = default_value;

			data = result;
		}
	}
};












#endif /* OP2A_READ_WRITE_HPP_ */
