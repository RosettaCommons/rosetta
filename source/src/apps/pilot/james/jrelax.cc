// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers

#include <protocols/jobdist/standard_mains.hh>
#include <basic/options/option.hh>

#include <devel/init.hh>

// C++ headers
#include <iostream>
#include <string>


// AUTO-REMOVED #include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <core/scoring/constraints/util.hh>
#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/relax/ClassicRelax.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {
	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::jobdist;
	using namespace protocols::moves;
	using namespace protocols::relax;
	using namespace core::scoring::constraints;

	ClassicRelax::register_options();
	FastRelax::register_options();
	register_options_universal_main();
	option.add_relevant( OptionKeys::in::file::fullatom );
	option.add_relevant( OptionKeys::relax::fast );
	devel::init(argc, argv);

	protocols::moves::MoverOP protocol = generate_relax_from_cmd();
	if ( option[ constraints::cst_fa_file ].user() ) {
		SequenceMoverOP seqmov = new SequenceMover;
		protocols::simple_moves::ConstraintSetMoverOP loadCsts( new protocols::simple_moves::ConstraintSetMover );
		loadCsts->constraint_file( get_cst_fa_file_option() );
		seqmov->add_mover( loadCsts );
		seqmov->add_mover( protocol );
		protocol = seqmov;
	}
	not_universal_main( *protocol );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
