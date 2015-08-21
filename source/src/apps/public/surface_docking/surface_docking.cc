// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   surface_docking.cc
/// @brief
/// @author Robin A Thottungal (raugust1@jhu.edu)
/// @author Michael Pacella (mpacella88@gmail.com)

// Unit header
#include <protocols/surface_docking/SurfaceDockingProtocol.hh>
#include <protocols/surface_docking/SurfaceDockingProtocol.fwd.hh>

// Project header
#include <protocols/jd2/JobDistributor.hh>
#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>
//#include <protocols/init/init.hh>

int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols::surface_docking;
		using namespace protocols::jd2;

		// initialize option system, random number generators, and all factory-registrators
		devel::init(argc, argv);
		//protocols::init(argc, argv);

		SurfaceDockingProtocolOP dp( new SurfaceDockingProtocol() );
		JobDistributor::get_instance()->go(dp);
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
