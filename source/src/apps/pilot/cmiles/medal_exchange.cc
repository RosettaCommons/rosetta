// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/cmiles/medal.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Utility headers
#include <devel/init.hh>

// Project headers
#include <protocols/medal/MedalMain.hh>
#include <protocols/viewer/viewers.hh>

//Auto Headers
#include <utility/vector1.hh>
int main(int argc, char* argv[]) {
	try {
		devel::init(argc, argv);
		protocols::viewer::viewer_main(
			protocols::medal::MedalExchange_main);
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
