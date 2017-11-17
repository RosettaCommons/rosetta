// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/cmiles/star_abinitio.cc
/// @author Christopher Miles (cmiles@uw.edu)

// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>

// Project headers
#include <protocols/star/StarAbinitio.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/viewer/viewers.hh>

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

	if ( !option[in::file::alignment].user() ) {
		utility_exit_with_message(ERROR_PREFIX + "-in:file:alignment");
	}

	if ( !option[in::file::template_pdb].user() ) {
		utility_exit_with_message(ERROR_PREFIX + "-in:file:template_pdb");
	}
}

void* star_main(void*) {
	using protocols::jd2::JobDistributor;
	using protocols::star::StarAbinitio;
	using protocols::moves::MoverOP;

	check_required();
	MoverOP mover( new StarAbinitio() );

	try {
		JobDistributor::get_instance()->go(mover);
	} catch (utility::excn::Exception& e) {
		std::cerr << "Exception: " << std::endl;
		e.show(std::cerr);
	}
	return 0;
}

int main(int argc, char* argv[]) {
	try {
		devel::init(argc, argv);
		protocols::viewer::viewer_main(star_main);
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
