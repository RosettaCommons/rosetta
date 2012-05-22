///////////////////////////////////////////////////////////////////////////////
//                                                                             
//  Copyright (C) 2010-2011  Artyom Beilis (Tonkikh) <artyomtnk@yahoo.com>     
//                                                                             
//  Distributed under:
//
//                   the Boost Software License, Version 1.0.
//              (See accompanying file LICENSE_1_0.txt or copy at 
//                     http://www.boost.org/LICENSE_1_0.txt)
//
//  or (at your opinion) under:
//
//                               The MIT License
//                 (See accompanying file MIT.txt or a copy at
//              http://www.opensource.org/licenses/mit-license.php)
//
///////////////////////////////////////////////////////////////////////////////
#ifndef CPPDB_NUMERIC_UTIL_H
#define CPPDB_NUMERIC_UTIL_H

#include <cppdb/errors.h>
#include <string>
#include <sstream>
#include <limits>
#include <iomanip>

namespace cppdb {

	///
	/// Small utility functions for backends, accepts - source string and stringstream with imbued std::locale
	/// it tries to case the value to T in best possible way.
	///
	/// For floating point string it casts it to the nearest ineger
	///
	template<typename T>
	T parse_number(std::string const &s,std::istringstream &ss)
	{
		ss.clear();
		ss.str(s);
		if(s.find_first_of(".eEdD")!=std::string::npos) {
			long double v;
			ss >> v;
			if(ss.fail() || !std::ws(ss).eof())
				throw bad_value_cast();
#ifdef __BORLANDC__
#pragma warn -8008 // condition always true/false
#pragma warn -8066 // unreachable code
#endif
			if(std::numeric_limits<T>::is_integer) {
				if(v > std::numeric_limits<T>::max() || v < std::numeric_limits<T>::min())
					throw bad_value_cast();
			}
#ifdef __BORLANDC__
#pragma warn .8008
#pragma warn .8066
#endif
			return static_cast<T>(v);
		}
		T v;
		ss >> v;
		if(ss.fail() || !std::ws(ss).eof()) 
			throw bad_value_cast();
		if(	std::numeric_limits<T>::is_integer 
			&& !std::numeric_limits<T>::is_signed 
			&& s.find('-') != std::string::npos 
			&& v!=0) 
		{
			throw bad_value_cast();
		}
		return v;
	}


}
#endif
