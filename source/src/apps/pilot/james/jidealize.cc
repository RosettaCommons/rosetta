// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Mike Tyka & Phil bradley
/// @brief


// libRosetta headers
#include <protocols/idealize.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/not_universal_main.hh>
#include <basic/options/option.hh>
#include <devel/init.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace core;

#include <basic/options/option_macros.hh>

#include <utility/excn/Exceptions.hh>

OPT_KEY( Real, atom_pair_constraint_weight )
OPT_KEY( Real, coordinate_constraint_weight )
OPT_KEY( Boolean, fast )

void
register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( atom_pair_constraint_weight, "atompair constraint weight", 0.0 );
	NEW_OPT( coordinate_constraint_weight, "coordinate constraint weight", 0.0 );
	NEW_OPT( fast, "fast protocol", false );
}

int
main( int argc, char * argv [] )
{
	try {

	using namespace protocols;
	using namespace protocols::jobdist;
	using namespace protocols::moves;
	using basic::options::option;
	using namespace basic::options::OptionKeys;


	register_options();
	devel::init(argc, argv);

	protocols::IdealizeMover idealizer;

	// set some options
	if ( option[ coordinate_constraint_weight ].user() ) {
		idealizer.coordinate_constraint_weight( option[ coordinate_constraint_weight ]() ) ;
	}
	if ( option[ atom_pair_constraint_weight ].user() ) {
		idealizer.atom_pair_constraint_weight( option[ atom_pair_constraint_weight ]() );
	}
	idealizer.fast( option[ fast ]() );

	not_universal_main( idealizer );

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

