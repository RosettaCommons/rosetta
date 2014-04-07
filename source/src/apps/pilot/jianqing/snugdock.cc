// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/pilot/jianqing/snugdock.cc
/// @brief
/// @author Jianqing Xu (xubest@gmail.com)
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)



#include <protocols/antibody/snugdock/SnugDockProtocol.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <utility/excn/Exceptions.hh>
//#include <utility/exit.hh>

#include <devel/init.hh>

int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols::antibody;
		using namespace protocols::jd2;

		SnugDockProtocol::register_options();
		protocols::jd2::register_options();

		// initialize core
		devel::init(argc, argv);

		SnugDockProtocolOP snugdock = new SnugDockProtocol;
		JobDistributor::get_instance()->go( snugdock );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
