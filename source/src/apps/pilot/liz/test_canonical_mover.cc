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
/// @author
/// @author

//

#include <protocols/canonical_sampling/CanonicalSamplingApplication.hh>
#include <protocols/canonical_sampling/CanonicalSamplingMover.hh>
#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>






int
main( int argc, char * argv [] ){
    try {

	protocols::canonical_sampling::CanonicalSamplingMover::register_options();

	devel::init(argc, argv);

	protocols::canonical_sampling::canonical_sampling_main();

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
}
