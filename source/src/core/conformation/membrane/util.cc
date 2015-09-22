// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// TODO: Get rid of this class

/// @file  core/conformation/membrane/util.fwd.hh
/// @brief  utility functions for membrane things
/// @author  JKLeman (julia.koehler1982@gmail.com)

// Unit headers
#include <core/conformation/membrane/util.hh>

// Project Headers
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>

// Package Headers
#include <core/types.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>

// Utility headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "core.conformation.membrane.util" );

namespace core {
namespace conformation {
namespace membrane {

// read spanfiles
utility::vector1< std::string > spanfile_names(){

	using namespace basic;
	using namespace basic::options;

	if ( ! option[OptionKeys::mp::setup::spanfiles].user() ) {
		utility_exit_with_message("Please provide spanfile(s)!");
	}

	// get filenames from Optionsystem
	utility::vector1< std::string > spanfiles = option[OptionKeys::mp::setup::spanfiles]();

	return spanfiles;

}// read spanfiles

// read single spanfile from options
std::string spanfile_name(){
	return spanfile_names()[1];
}

} // membrane
} // conformation
} // core
