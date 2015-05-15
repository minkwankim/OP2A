/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 11, 2015
 *      			Author: Minkwan Kim
 *
 * Version.hpp
 * 			-  
 *  
 */
#ifndef VERSION_HPP_
#define VERSION_HPP_




#include "Common/include/Standard_headers.hpp"
#include "Common/include/Common.hpp"


namespace OP2A {
namespace Common {
	/*
	 * Version information.
	 */

	class Common_API Version
	{
		public:
		  explicit Version (int ver_main, int ver_sub, int date_year, int date_month, int date_day, const char * type_info);

		  void info() const;

		private:
			int primary;       		/* Primary version number   (i.e. 2 for V2.0) */
			int secondary;    	 	/* Secondary version number (i.e. 0 for V2.0) */
			int year;          		/* Most recent date of update */
			int month;
			int date;

			std::string type;     	/* "Type" of version */
	};
}
}


#endif /* VERSION_HPP_ */
