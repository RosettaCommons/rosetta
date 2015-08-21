// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/star/StarAbinitioMain.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/star/StarAbinitioMain.hh>

// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/exit.hh>
#include <utility/excn/EXCN_Base.hh>

// Project headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>

// Package headers
#include <protocols/star/StarAbinitio.hh>

namespace protocols {
namespace star {

const std::string ERROR_PREFIX = "Failed to specify required option ";

void check_required() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( !option[in::file::fasta].user() ) {
		utility_exit_with_message(ERROR_PREFIX + "-in:file:fasta");
	}

	if ( !option[in::file::frag3].user() ) {
		utility_exit_with_message(ERROR_PREFIX + "-in:file:frag3");
	}

	if ( !option[in::file::frag9].user() ) {
		utility_exit_with_message(ERROR_PREFIX + "-in:file:frag9");
	}

	if ( !option[cm::aln_format].user() ) {
		utility_exit_with_message(ERROR_PREFIX + "-cm:aln_format");
	}

	if ( !option[in::file::alignment].user() ) {
		utility_exit_with_message(ERROR_PREFIX + "-in:file:alignment");
	}

	if ( !option[in::file::template_pdb].user() ) {
		utility_exit_with_message(ERROR_PREFIX + "-in:file:template_pdb");
	}
}

void* graphics_main(protocols::moves::MoverOP mover) {
	using protocols::jd2::JobDistributor;

	try {
		JobDistributor::get_instance()->go(mover);
	} catch (utility::excn::EXCN_Base& e) {
		std::cerr << "Exception: " << std::endl;
		e.show(std::cerr);
	}
	return 0;
}

void* StarAbinitio_main(void*) {
	check_required();
	return graphics_main(protocols::moves::MoverOP( new StarAbinitio() ));
}

}  // namespace star
}  // namespace protocols
