// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief   Setup a membrane protein by adding the membrane virtual
/// @author  JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <protocols/membrane/AddMembraneMover.hh>

// Package Headers
#include <apps/benchmark/performance/init_util.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <basic/Tracer.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using basic::Error;
using basic::Warning;

using namespace core;
using namespace core::pose;
using namespace core::conformation;
using namespace core::conformation::membrane;

static basic::Tracer TR( "apps.pilot.jkleman.add_membrane" );

int
main( int argc, char * argv [] ) {
	try {

		using namespace protocols::jd2;
		using namespace protocols::membrane;

		devel::init(argc, argv);

		// create two movers and concatenate them in a sequence mover
		AddMembraneMoverOP addmem( new AddMembraneMover() );

		// call jobdistributor on sequence mover
		JobDistributor::get_instance()->go( addmem );
	}
catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}

}
