// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file combine_silent.cc
/// @brief simple application for combining a number of silent-files into a
/// single silent-file.
/// @author James Thompson

// libRosetta headers

#include <devel/init.hh>

#include <basic/options/option.hh>

#include <protocols/moves/NullMover.hh>
#include <protocols/jobdist/not_universal_main.hh>

// C++ headers
#include <string>


// option key includes

#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


int
main( int argc, char* argv [] ) {
	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		// define relevant options
		OPT( in::path::database);
		OPT( in::file::s );
		OPT( in::file::l );
		OPT( in::file::silent );
		OPT( in::file::tags );
		OPT( in::file::silent_struct_type );
		OPT( in::file::silent_renumber );
		OPT( in::file::residue_type_set );
		OPT( out::file::silent );
		OPT( out::file::silent_struct_type );

		// options, random initialization
		devel::init( argc, argv );

		std::string usage("");
		usage += "\n\nusage:  combine_silent [options] -in::file::silent <silent_files> or -in::file::s <pdb> or -in::file::l <list of pdb files> \n";
		usage += "\tTo see a list of other valid options, use the option -help.\n";

		if ( !option[ in::file::silent ].user() && !option[ in::file::s ].user() && !option[ in::file::l ].user() ) {
			std::cerr << usage << std::endl;
			std::exit(1);
		}

		protocols::moves::NullMover mover;
		protocols::jobdist::not_universal_main( mover );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // main
