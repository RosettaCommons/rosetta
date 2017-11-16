// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/smlewis/jd2memtest
/// @brief empty app to check jd2 memory usage against many files
/// @author Steven Lewis


// Unit Headers


// Project Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/NullMover.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/Tracer.hh>


int main( int argc, char* argv[] )
{
	try {

		devel::init(argc, argv);

		protocols::jd2::JobDistributor::get_instance()->go(new protocols::moves::NullMover);

		T("jd2memtest") << "************************d**o**n**e**************************************" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
