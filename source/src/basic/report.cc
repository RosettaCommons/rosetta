// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/report.cc
/// @brief  Report and Reporter
/// @author Sergey Lyskov

#include <basic/report.hh>

#include <utility/json_spirit/json_spirit_writer.h>

#include <fstream>

using std::string;

namespace basic {

Report::Report(string const & file_name)
{
	file_name_ = file_name;
}

Report::~Report()
{
	write();
}


void Report::write()
{
	std::ofstream ft( (file_name_+".txt").c_str() ); ft << text_;
	std::ofstream fj( (file_name_+".json").c_str() ); fj << utility::json_spirit::write(data_, utility::json_spirit::pretty_print);
}



} // namespace basic
