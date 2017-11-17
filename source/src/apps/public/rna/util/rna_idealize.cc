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

#include <protocols/rna/movers/RNAIdealizeMover.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/option.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/viewer/viewers.hh>
#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::rna::movers;

OPT_KEY( Integer, iterations )
OPT_KEY( Boolean, noise )
OPT_KEY( Boolean, final_minimize )
OPT_KEY( Real, ang_significance_threshold )
OPT_KEY( Boolean, handle_suites )


void*
my_main ( void* ) {

	RNAIdealizeMoverOP rim( new RNAIdealizeMover( option[ iterations ].value(), option[ noise ].value(), option[ final_minimize ].value(), option[ ang_significance_threshold ].value(), option[ handle_suites ].value()  ) );

	protocols::jd2::JobDistributor::get_instance()->go( rim );

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}

///////////////////////////////////////////////////////////////////////////////
int
main ( int argc, char * argv [] ) {
	try {
		NEW_OPT( iterations, "number of iterations across which to spread idealization", 100 );
		NEW_OPT( noise, "add noise to the initial coordinates, providing a useful reason to run this more than once", false );
		NEW_OPT( final_minimize, "do a final, unrestrained minimize", false );
		NEW_OPT( ang_significance_threshold, "Size of angle deviation to correct (degrees)", 5 );
		NEW_OPT( handle_suites, "repair suite outliers", false );

		devel::init( argc, argv );

		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
