/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * ProcessInfo.cpp
 * 			-  
 *  
 */

#include <sstream>
#include "Common/include/ProcessInfo.hpp"


using namespace std;

namespace OP2A{
namespace Common{


ProcessInfo::ProcessInfo()
{

}

ProcessInfo::~ProcessInfo()
{

}

string ProcessInfo::memoryUsage() const
{
	const OPdouble bytes = memoryUsageBytes();

	ostringstream out;
	if (bytes/1024 <= 1)				out << bytes << "B";
	else if (bytes/1024/1024 <= 1)		out << bytes << "kB";
	else if (bytes/1024/1024/2014 <= 1)	out << bytes << "MB";
	else 								out << bytes << "GB";

	return out.str();
}


} // end namespace Common
} // end namespace OP2A


