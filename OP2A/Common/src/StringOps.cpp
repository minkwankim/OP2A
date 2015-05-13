/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 13, 2015
 *      			Author: Minkwan Kim
 *
 * StringOps.cpp
 * 			-  
 *  
 */

#include <cctype>
#include "Common/include/StringOps.hpp"



using namespace std;

namespace OP2A {
namespace Common {



void StringOps::join (const std::string& delim, std::string parts[], unsigned int nparts, std::string& out)
{
	out = parts[0];
	for (size_t i = 1; i < nparts; ++i)
	{
		out += delim;
		out += parts[i];
	}
}


void StringOps::toLower (std::string& out)
{
	const string::iterator end = out.end();
	for (string::iterator itr = out.begin(); itr != end; ++itr)
	{
		*itr = tolower(*itr);
	}
}


void StringOps::toUpper (std::string& out)
{
	const string::iterator end = out.end();
	for (string::iterator itr = out.begin(); itr != end; ++itr)
	{
		*itr = toupper(*itr);
	}
}


void StringOps::subst (const std::string& lhs, const std::string& rhs, std::string& out)
{
	std::string::size_type i = 0;
	while (i < out.length())
	{
		std::string substr (out, i);
		std::string::size_type idx = out.find(lhs, i);
		if (idx == std::string::npos)
		{
			break;
		}
		else
		{
			out.replace(idx, lhs.length(), rhs);
			i = idx + rhs.length();
			if (!rhs.length()) ++i;
		}
	}
}


void StringOps::trimFront (std::string& out)
{
	size_t startpos = out.find_first_not_of(" \t");
	if( string::npos != startpos )	out = out.substr( startpos );
}


void StringOps::trimRear (std::string& out)
{
	size_t endpos = out.find_last_not_of(" \t");
	if( string::npos != endpos ) out = out.substr( 0, endpos+1 );
}


void StringOps::trim (std::string& out)
{
	size_t startpos = out.find_first_not_of(" \t");
	size_t endpos = out.find_last_not_of(" \t");

	if(( string::npos == startpos ) || ( string::npos == endpos))	out = "";
	else															out = out.substr( startpos, endpos-startpos+1 );
}


void StringOps::trim2 (std::string& out)
{
	size_t startpos = out.find_first_not_of(" \t\n");
	size_t endpos = out.find_last_not_of(" \t\n");

	if(( string::npos == startpos ) || ( string::npos == endpos))	out = "";
	else															out = out.substr( startpos, endpos-startpos+1);
}


std::vector<std::string> StringOps::getWords (const std::string& in)
{
	std::string s = in;            // copy
	vector<std::string> words;     // return vector

	bool inWord = false;           // whether we're in a word
	std::string word;              // current word
	for (size_t i = 0, len = s.length(); i < len; ++i)
	{
		OPchar ch = s[i];
		if (inWord)
		{
			if (isspace(ch))
			{
				words.push_back(word);
				word = "";
				inWord = false;
			}
			else
			{
				word += ch;
			}
		}
		else
		{
			if (!isspace(ch))
			{
				word += ch;
				inWord = true;
			}
		}
	}

	if (inWord)	words.push_back(word);
	return words;
}


vector<std::string> StringOps::getWords  (const std::string& in, const OPchar sep )
{
	std::string s = in;                // copy
	vector<std::string> words;         // return vector

	bool inWord = false;        // whether we're in a word
	std::string word;                // current word
	for (size_t i = 0, len = s.length(); i < len; ++i)
	{
		OPchar ch = s[i];
		if (inWord)
		{
			if (ch == sep)
			{
				words.push_back(word);
				word = "";
				inWord = false;
			}
			else
			{
				word += ch;
			}
		}
		else if (ch != sep)
		{
			word += ch;
			inWord = true;
		}
	}

	if (inWord)	words.push_back(word);
	return words;
}


bool StringOps::startsWith (const std::string& in, const std::string& str)
{
	return in.length() >= str.length() && in.substr(0, str.length()) == str;
}


bool StringOps::endsWith (const std::string& in, const std::string& str)
{
	return in.length() >= str.length() && in.substr(in.length() - str.length()) == str;
}

std::string StringOps::convertYesNo(const int i)
{
	std::string ss;

	if (i == 1) ss = "YES";
	else if (i == 0) ss = "NO";

	return ss;
}


} // namespace Common
} // namespace COOLFluiD
