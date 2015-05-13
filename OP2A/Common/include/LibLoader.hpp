/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 12, 2015
 *      			Author: Minkwan Kim
 *
 * LibLoader.hpp
 * 			-  
 *  
 */
#ifndef LIBLOADER_HPP_
#define LIBLOADER_HPP_

#include <boost/filesystem/path.hpp>

#include "Common/include/Common.hpp"
#include "Common/include/NonCopyable.hpp"
#include "Common/include/Exception.hpp"
#include "Common/include/Exception_LibLoader.hpp"


namespace OP2A{
namespace Common{

class Common_API LibLoader:public Common::NonCopyable<LibLoader>
{
public:
	LibLoader(){};
	virtual ~LibLoader(){};

	virtual void load_library(const std::string& lib) = 0;
	virtual void set_search_paths(std::vector< boost::filesystem::path >& paths) = 0;


	static std::string getClassName() { return "LibLoader"; }
};

}
}



#endif /* LIBLOADER_HPP_ */
