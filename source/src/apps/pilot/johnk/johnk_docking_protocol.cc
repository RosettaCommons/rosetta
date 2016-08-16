// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @ docking_protocol.cc
/// @ author James Thompson, Jeff Gray


// Rosetta headers
#include <protocols/docking/DockingProtocol.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>

#include <devel/init.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>
//#include <utility/tools/make_vector1.hh>

// option key includes
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <core/scoring/ScoreFunctionFactory.hh>


// JQX added below two lines for forcing the code
// to use the constant seed in production mode
// Sergey believes that in production run,
// the -constant_seed is not working properly in Jump.cc file
// Please see the details in the Jump.cc file
//#include <core/init.hh>     //JQX
//#include <numeric/random/random.hh> //JQX


int
main( int argc, char * argv [] )
{
	using namespace basic::options;
	using namespace protocols::docking;
	using namespace protocols::jd2;

	DockingProtocol::register_options();
	protocols::jd2::register_options();
	// initialize core
	devel::init(argc, argv);
//	core::init_random_generators(3,numeric::random::_RND_TestRun_, "mt19937"); //JQX from Sergery

	core::scoring::ScoreFunctionOP highres_scorefxn = core::scoring::get_score_function();

//	DockingProtocolOP dp = new DockingProtocol( utility::tools::make_vector1<core::SSize>(1), false, false, true, NULL, highres_scorefxn );
	DockingProtocolOP dp = new DockingProtocol( 1, false, true, true, NULL, highres_scorefxn );

	if ( option[ OptionKeys::docking::multibody ].user() ) {
			utility::vector1< core::Size > const movable_jumps(
				option[ OptionKeys::docking::multibody ]()
			);
			dp->set_movable_jumps( movable_jumps );
	}

	JobDistributor::get_instance()->go(dp);
}

