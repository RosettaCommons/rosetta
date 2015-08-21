// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @ docking_Prepackprotocol.cc
/// @ author Robin A Thottungal


// Unit header
#include <protocols/docking/DockingPrepackProtocol.hh>
#include <protocols/docking/DockingEnsemblePrepackProtocol.hh>
#include <protocols/docking/DockingHighRes.fwd.hh>

// Project header
#include <protocols/jd2/JobDistributor.hh>
#include <devel/init.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


int
main( int argc, char * argv [] )
{
	try {

		using namespace basic::options;
		using namespace protocols::docking;
		using namespace protocols::jd2;

		// initialize core
		devel::init(argc, argv);
		DockingHighResOP dp;

		if ( option[ OptionKeys::docking::ensemble1 ].user() || option[ OptionKeys::docking::ensemble2 ].user() ) {
			dp = DockingHighResOP( new DockingEnsemblePrepackProtocol() );
		} else dp = DockingHighResOP( new DockingPrepackProtocol() );

		JobDistributor::get_instance()->go(dp);
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

