/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: Jun 22, 2015
 *      			Author: Minkwan Kim
 *
 * DataReadFromFile.cpp
 * 			-  
 *  
 */

#include "Common/include/DataReadFromFile.hpp"


using namespace std;


namespace OP2A {
namespace Common {


// ::F-01
// Read date from string data
std::string DataRead::get_data_from_string_string(std::string line, std::string TEXT)
{
	int num;
	std::string result;

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
		result = line;
	}
	else
	{
		result = "";
	}

	return (result);
}



// ::F-02
// Remove front space in string data
string DataRead::remove_space_front(string line)
{
	int i, num;
	string result;

	num = line.size();
	for (i = 0; i <= num-1; i++)
	{
		if(line.compare(0, 1," ") == 0) line.erase(0, 1);
	}

	result = line;
	return result;
}



// ::F-03
// Remove end space in string data
string DataRead::remove_space_end(string line)
{
	int i, num, num_temp;
	string result;

	num = line.size();
	for (i = 0; i <= num-1; i++)
	{
		num_temp = line.size();
		if(line.compare(num_temp-1, 1," ") == 0) line.erase(num_temp-1, 1);
	}
	result = line;
	return result;
}


// ::F-04
// Remove comments
void DataRead::remove_comments(string& line, string comment_s, string comment_e)
{
	string result_text;

	int i;
	int size1, size2, size3;
	int index_s, index_e, index_size;
	size1 = line.size();
	size2 = comment_s.size();
	size3 = comment_e.size();


	index_s = -1;
	index_e = -1;
	index_size = 0;
	for (i = 0; i <= size1-1; i++)
	{
		if (line.compare(i, size2, comment_s) == 0)
		{
			index_s = i;
		}

		if (line.compare(i, size3, comment_e) == 0)
		{
			index_e = i;
		}
	}

	if (index_s != -1)
	{
		if (index_e != -1) index_size = index_e - index_s + size3;
		else index_size = size1 - index_s;

		line.erase(index_s, index_size);
	}
	result_text = line;

	// return result_text;
}

void DataRead::remove_comments(string& line, string comment_s)
{
	string result_text;

	int i;
	int size1, size2;
	int index_s, index_size;
	size1 = line.size();
	size2 = comment_s.size();


	index_s = -1;
	index_size = 0;
	for (i = 0; i <= size1-1; i++)
	{
		if (line.compare(i, size2, comment_s) == 0)
		{
			index_s = i;
		}
	}

	if (index_s != -1)
	{
		index_size = size1 - index_s;
		line.erase(index_s, index_size);
	}
	result_text = line;

	//return result_text;
}



// ::F-05
// Read line
void DataRead::read_line(string& line)
{
	/////////////////////////////////////////////////////
	// REMOVE COMMENTS: SUPPORT C/C++ AND MATLAB STYLE //
	/////////////////////////////////////////////////////
	// C/C++ STYLE
	remove_comments(line, "/*", "*/");
	remove_comments(line, "//");

	// MATLAB STYLE
	remove_comments(line, "/%", "%/");
	remove_comments(line, "%");
}








}
}
