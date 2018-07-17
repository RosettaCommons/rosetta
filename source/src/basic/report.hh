// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/report.hh
/// @brief  Report and Reporter
/// @author Sergey Lyskov


#ifndef INCLUDED_basic_Report_hh
#define INCLUDED_basic_Report_hh

#include <basic/report.fwd.hh>

#include <utility/json_spirit/json_spirit_value.h>
#include <utility/json_spirit/json_spirit_writer_options.h>

#include <string>
#include <sstream>

namespace basic {

/// @brief Data class to hold Report information. Use Reporter class for write access
class Report
{
	std::string file_name_;

	/// @brief store human readable text info
	std::string text_;

	/// @brief store json data
	utility::json_spirit::Object data_;

public:
	Report(std::string const & file_name);
	~Report();

	void write();

	template <typename T>
	Report& operator<<(T const &v) { std::ostringstream s; s << v;  text_+=s.str();  return *this; }

	template <typename T>
	void set(std::string const &key, T const &value) { data_.push_back( utility::json_spirit::Pair(key, value) ); }

};


/// @brief Proxy class to access Reporter class though OP. Allow null OP on init (all reporter function in this case do nothing)
class Reporter
{
	ReportOP report_;

public:
	Reporter(ReportOP report=basic::ReportOP()) { report_ =  report; }

	virtual ~Reporter() {};

	template <typename T>
	Reporter& operator<<(T const &v) { if ( report_ ) *report_ << v;
		return *this; }

	template <typename T>
	void set(std::string const &key, T const &value) { if ( report_ ) report_->set(key, value); }
};


} // namespace basic

#endif // INCLUDED_basic_Report_hh
