// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Mike Tyka & Phil Bradley
/// @brief


// libRosetta headers
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <basic/options/option.hh>

#include <devel/init.hh>



// C++ headers
//#include <cstdlib>
// AUTO-REMOVED #include <fstream>
#include <iostream>
#include <string>

using namespace core;

#include <basic/options/option_macros.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

#include <protocols/jobdist/Jobs.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


OPT_KEY( Real, atom_pair_constraint_weight )
OPT_KEY( Real, coordinate_constraint_weight )
OPT_KEY( Boolean, fast )
OPT_KEY( Boolean, chainbreaks )

void
register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( atom_pair_constraint_weight, "atompair constraint weight", 0.0 );
	NEW_OPT( coordinate_constraint_weight, "coordinate constraint weight", 0.0 );
	NEW_OPT( fast, "fast protocol", false );
	NEW_OPT( chainbreaks, "keep chainbreaks", false );
}

int
main( int argc, char * argv [] )
{
	using namespace protocols;
	using namespace protocols::moves;
	using basic::options::option;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	register_options();
	devel::init(argc, argv);

	// set some idealize mover
	protocols::idealize::IdealizeMoverOP idealizer = new protocols::idealize::IdealizeMover;
	if ( option[ coordinate_constraint_weight ].user() ) {
		idealizer->coordinate_constraint_weight( option[ coordinate_constraint_weight ]() ) ;
	}
	if ( option[ atom_pair_constraint_weight ].user() ) {
		idealizer->atom_pair_constraint_weight( option[ atom_pair_constraint_weight ]() );
	}
	idealizer->fast( option[ fast ]() );
	idealizer->chainbreaks( option[ chainbreaks ]() );

	MoverOP mover (idealizer);

	// optionally set pose for symmetry
	if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
		protocols::moves::SequenceMoverOP seqmov = new protocols::moves::SequenceMover;
		seqmov->add_mover( new protocols::simple_moves::symmetry::SetupForSymmetryMover );
		seqmov->add_mover( mover );
		mover = seqmov;
	}

	protocols::jobdist::main_plain_mover( *mover );

	return 0;
}

