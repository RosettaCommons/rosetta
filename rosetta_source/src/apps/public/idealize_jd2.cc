// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Mike Tyka
/// @author Phil bradley
/// @author Christopher Miles (cmiles@uw.edu)
/// @brief
#include <protocols/idealize/idealize.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/idealize/IdealizeMover.fwd.hh>
#include <protocols/moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/options/option.hh>
#include <devel/init.hh>

// Utility
#include <basic/options/option_macros.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

// C++ headers
#include <iostream>

OPT_KEY( Real, atom_pair_constraint_weight )
OPT_KEY( Real, coordinate_constraint_weight )
OPT_KEY( Boolean, fast )

void register_options() {
  NEW_OPT( atom_pair_constraint_weight, "atompair constraint weight", 0.0 );
  NEW_OPT( coordinate_constraint_weight, "coordinate constraint weight", 0.0 );
  NEW_OPT( fast, "fast protocol", false );
}

int main( int argc, char * argv [] ) {
	using namespace protocols::idealize;
	using namespace protocols::jd2;
  using namespace protocols::moves;
  using namespace protocols::moves::symmetry;
	using namespace basic::options;
  using namespace basic::options::OptionKeys;

	// initialization
  register_options();
  devel::init(argc, argv);

	// configure the idealize mover
	IdealizeMoverOP idealizer = new IdealizeMover();

	if ( option[ coordinate_constraint_weight ].user() ) {
		idealizer->coordinate_constraint_weight( option[ coordinate_constraint_weight ]() ) ;
	}

	if ( option[ atom_pair_constraint_weight ].user() ) {
		idealizer->atom_pair_constraint_weight( option[ atom_pair_constraint_weight ]() );
	}
	idealizer->fast( option[ fast ]() );

	MoverOP mover (idealizer);

	// optionally set pose for symmetry
	if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
		SequenceMoverOP seqmov = new SequenceMover();
		seqmov->add_mover( new SetupForSymmetryMover() );
		seqmov->add_mover( mover );
		mover = seqmov;
	}

	// start the job
  try {
    JobDistributor::get_instance()->go( mover );
  } catch ( utility::excn::EXCN_Base& excn ) {
    std::cerr << "Exception: " << std::endl;
    excn.show( std::cerr );
  }
}
