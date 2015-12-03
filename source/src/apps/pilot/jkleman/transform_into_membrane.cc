// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief   Transform a protein into the membrane
/// @author  JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/moves/MoverContainer.hh>

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

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.jkleman.transform_into_membrane" );

int
main( int argc, char * argv [] ) {
	try {

		using namespace protocols::jd2;
		using namespace protocols::membrane;
		using namespace protocols::moves;

		devel::init(argc, argv);

		// create two movers and concatenate them in a sequence mover
		AddMembraneMoverOP addmem( new AddMembraneMover() );
		TransformIntoMembraneMoverOP transform( new TransformIntoMembraneMover() );
		SequenceMoverOP seq( new SequenceMover( addmem, transform ) );

		// call jobdistributor on sequence mover
		JobDistributor::get_instance()->go( seq );
	}
catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}

}
