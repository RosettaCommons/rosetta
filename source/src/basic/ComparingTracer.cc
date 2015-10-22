// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/basic/ComparingTracer.hh
/// @brief  version of Tracer calss for to tracking down the source of an instability by compare
///         output with file.
/// @author Sergey Lyskov


#include <basic/ComparingTracer.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <iostream>
#include <cassert>
#include <cstddef>
#include <iosfwd>

namespace basic {

ComparingTracer::ComparingTracer(std::string const & file_name)
{
	using namespace basic;
	using namespace basic::options::OptionKeys;

	dry_run_ = basic::options::option[ out::dry_run ]();

	if ( dry_run_ ) {
		file_.open(file_name.c_str(), std::fstream::out);
	} else {
		file_.open(file_name.c_str(), std::fstream::in);
	}
}

ComparingTracer::~ComparingTracer()
{
	// We want to insure that last symbol in output is a new line,
	// and that inner buffer is flushed.
	(*this) << std::endl;
	file_.close();
}


void ComparingTracer::t_flush(std::string const & s)
{
	if ( dry_run_ ) {
		file_ << s;
	} else {
		for ( size_t i=0; i<s.size(); i++ ) {
			if ( file_.get() != s[i] ) {
				std::cerr << "ComparingTracer:: Output is differnet from the original: " << s << std::endl;
				assert(false);
			}
		}
	}
}


} // namespace basic

