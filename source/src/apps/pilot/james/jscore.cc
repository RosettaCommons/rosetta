// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file jscore.cc
/// @brief

#include <core/scoring/constraints/util.hh>

#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/NullMover.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>

#include <devel/init.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <utility/excn/Exceptions.hh>

//static basic::Tracer TR("apps.jscore");

std::string get_env_var( std::string const & key ) {
	char * val;
	val = getenv( key.c_str() );
	std::string retval = "";
	if (val != NULL) {
		retval = val;
	}
	return retval;
}

int
main( int argc, char * argv [] ) {
	try {

	using namespace protocols::jobdist;
	using namespace protocols::moves;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//option[ in::file::rescore ].def( true );
	//option[ in::file::residue_type_set ].def( "fa_standard" );
	//option[ in::file::silent_struct_type ].def( "binary" );
	//option[ score::weights ].def( "score" );
	//option[ score::patch ].def( "score12" );
	//option[ out::file::silent_struct_type ].def( "score" );
	//option[ out::file::silent ].def( "scores.sc" );

	// initialize core
	devel::init(argc, argv);

	std::cout << "currently in file " << __FILE__ << std::endl;
	std::cout << "currently at line " << __LINE__ << std::endl;

	MoverOP mover( new NullMover );
	if ( option[ OptionKeys::constraints::cst_weight ].user() ) {
		protocols::simple_moves::ConstraintSetMoverOP loadcsts( new protocols::simple_moves::ConstraintSetMover );
		loadcsts->constraint_file( core::scoring::constraints::get_cst_file_option() );
		mover = loadcsts;
	}
	not_universal_main( *mover );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // main
